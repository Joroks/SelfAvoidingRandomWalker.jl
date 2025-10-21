struct Particle
    pos::Vector3
    shape::Matrix3x3
    pvalue::Vector3

    Particle(position::AbstractVector, radius::Real) = new(position, radius*I(3), [2, 2, 2])
    Particle(position::AbstractVector, shape::AbstractArray) = new(position, shape, [2, 2, 2])
    Particle(position::AbstractVector, shape::AbstractArray, pvalue::Real) = new(position, shape, ones(3)*pvalue)
    Particle(position::AbstractVector, radius::Number, pvalue::Real) = new(position, radius*I(3), ones(3)*pvalue)
    Particle(position::AbstractVector, shape::AbstractArray, pvalue::AbstractVector) = new(position, shape, pvalue)
end

mutable struct SARW{CG}
    box::Matrix3x3
    cell::Matrix3x3
    num_cells::NTuple{3, Int}
    dmin::Float64
    max_tries::Int

    chain_generator::CG

    particles::Memory{Particle}
    atoms::Memory{Vector3}
    tries::Memory{Int}
    chains::Memory{Int}
    cells::CountingBins{3}
    
    chain_counter::Int

    isdone::Bool
    highest_count::Int
    all_tries::Int

    function SARW(box, dmin, chains, max_tries, chain_generator::CG, particles=Particle[]) where CG
        self = new{CG}()

        self.box = Matrix3x3(box)
        self.num_cells, self.cell = calculate_cell(self.box, dmin)
        self.dmin = dmin
        self.max_tries = max_tries
        
        self.chain_generator = chain_generator

        natoms = sum(chains)
        self.particles = Memory{Particle}(undef, length(particles))
        self.atoms = Memory{Vector3}(undef, natoms)
        self.tries = Memory{Int}(undef, natoms)
        self.chains = Memory{Int}(undef, length(chains))
        self.cells = CountingBins(natoms, self.num_cells)

        self.particles .= particles
        self.tries .= 0
        cumsum!(self.chains, chains)
        self.chain_counter = 1
        self.isdone = false
        self.highest_count = 0
        self.all_tries = 0
        return self
    end
end

function inside_particle(self::SARW, particle::Particle, point::Vector3)
    diff = minimum_image(self, particle.pos - point)
    diff = particle.shape \ diff
    mapreduce(+, diff, particle.pvalue) do d, p
        abs(d)^p
    end <= 1
end


function calculate_cell(box::SMatrix, dmin)
    M = map(eachrow(inv(box))) do a
        n = inv(norm(a)*dmin)
        floor(Int, n)
    end
    C = box / Diagonal(M)
    return Tuple(M), C
end

function inscribed_circles(box)
    map(eachrow(inv(box))) do a
        1/norm(a)
    end
end

current_atom(self::SARW) = count(self.cells)+1
function current_segment(self::SARW, n=1)
    chain_start = self.chain_counter == 1 ? 1 : self.chains[self.chain_counter-1]+1
    atom = count(self.cells)
    max(chain_start, atom+1-n):atom
end

function continue_chain!(self::SARW)
    unit_vector(α, β) = SA[-cos(α), -sin(α)*sin(β), sin(α)*cos(β)]
    random_unit_vector() = unit_vector(acos(2rand()-1), 2pi*rand())

    atom = @view self.atoms[current_atom(self)]
    segment = @view self.atoms[current_segment(self, 3)]
    try_count = self.tries[current_atom(self)] += 1

    length(segment) == 0 && return atom[] = self.box*rand(Vector3)

    try_factor = (try_count - rand()) / self.max_tries
    (bond, angle, torsion) = self.chain_generator(try_factor)

    length(segment) == 1 && return atom[] = segment[end] + bond*random_unit_vector()

    r1 = normalize(segment[end] - segment[end-1])
    r2 = length(segment) == 2 ? random_unit_vector() : segment[2] - segment[1]
    n = normalize(r1 × r2)

    return atom[] = segment[end] + bond*[r1 n r1×n]*unit_vector(angle, torsion)
end

function backtrack!(self::SARW)
    while self.tries[current_atom(self)] >= self.max_tries
        isempty(current_segment(self)) && return
        self.tries[current_atom(self)] = 0
        decrement!(self.cells)
    end
end

function check_atom_position!(self::SARW)
    atom = self.atoms[current_atom(self)]
    prev = current_segment(self, 1)
    I = cell_index(self, atom)
    for offset in CartesianIndices((-1:1, -1:1, -1:1))
        J = index_mod(self, I+offset)
        for neighbor_index in self.cells[J]
            neighbor_index in prev && continue
            neighbor = self.atoms[neighbor_index]
            if periodic_distance(self, atom, neighbor) < self.dmin
                backtrack!(self)
                return
            end
        end
    end

    for particle in self.particles
        inside_particle(self, particle, atom) || continue
        backtrack!(self)
        return
    end

    if current_atom(self) >= self.chains[self.chain_counter]
        self.chain_counter += 1
    end

    increment!(self.cells[I])
end

function _run_SARW!(self::SARW)
    while !self.isdone
        continue_chain!(self)
        check_atom_position!(self)

        self.highest_count = max(self.highest_count, count(self.cells))
        self.all_tries += 1
        self.isdone |= self.chain_counter > length(self.chains)
    end
end

function run_SARW!(self::SARW)
    p = Progress(last(self.chains);
        desc = "Running SARW...",
        showspeed=true
    )

    function showvalues()
        speed = (time() - p.core.tinit)/self.all_tries
        SA[
            ("atoms placed", self.highest_count),
            ("number of candidates", self.all_tries),
            ("average number of candidates", self.all_tries/self.highest_count),
            ("time per candidate", ProgressMeter.speedstring(speed))
        ]
    end

    try @sync begin
        worker = Threads.@spawn _run_SARW!(self)

        while !istaskdone(worker)
            update!(p, self.highest_count, showvalues=showvalues())
            sleep(p.dt)
        end
    end finally
        self.isdone = true
    end

    finish!(p; showvalues=showvalues())
end

function index_mod(self::SARW, I::CartesianIndex)
    J = mod.(I.I, axes(self.cells))
    CartesianIndex(J)
end

function cell_index(self::SARW, atom::Vector3)
    I = floor.(Int, self.cell\atom) .+ 1
    index_mod(self, CartesianIndex(I...))
end

function minimum_image(self::SARW, point)
    point - self.box * round.(self.box \ point)
end

function periodic_distance(self::SARW, a, b)
    diff = minimum_image(self, a - b)
    norm(diff)
end

function atom_positions(self::SARW)
    current_start = 1

    map(self.chains) do l
        range = current_start:l
        current_start = l+1
        self.atoms[range]
    end
end
