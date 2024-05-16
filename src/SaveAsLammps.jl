using Printf
using ProgressMeter

"""
    saveAsLammps(boxSize, chains, filePath, headerTitle)

Save the result of the SARW algorithm to a lammps data file.

There are certain limitations on `boxSize` in order to generate a valid lammps file. It has to have the following structure:
```
    xhi  xy   xz
    0    yhi  yz
    0    0    zhi
```
with
```
    xhi, yhi, zhi > 0

    -xhi/2 ≤ xy ≤ xhi/2
    -xhi/2 ≤ xz ≤ yhi/2
    -xhi/2 ≤ yz ≤ yhi/2
```

for more information see: https://docs.lammps.org/Howto_triclinic.html

Note that any box size can be transformed into the valid structure using lattice reduction and QR-decomposition while retaining eqivalent periodicity (although prossibly rotated).

"""
function saveAsLammps(boxSize, chains, filePath, headerTitle)
    atoms = Tuple{Int64, SVector{3, Float64}}[]
    bonds = Tuple{Int64, Int64}[]
    angles = Tuple{Int64, Int64, Int64}[]
    

    for (chainIndex, chain) in enumerate(chains)
        for (atomIndex, atom) in enumerate(chain)
            Base.push!(atoms, (chainIndex, atom))
            numAtoms = length(atoms)

            if atomIndex >= 2
                Base.push!(bonds, (numAtoms-1, numAtoms))
            end

            if atomIndex >= 3
                Base.push!(angles, (numAtoms-2, numAtoms-1, numAtoms))
            end
        end
    end

    open(filePath, "w") do file
        desc = "Writing to LAMMPS file:"
        barlen = ProgressMeter.tty_width(desc, stderr, true) - 1 # workaround for incorrect length calculation
        p = Progress(length(atoms) + length(bonds)+ length(angles), desc, barlen=barlen)

        println(file, headerTitle)
        println(file)

        @printf file "%d atoms\n" length(atoms)
        @printf file "%d bonds\n" length(bonds)
        @printf file "%d angles\n\n" length(angles)

        println(file, "1 atom types")
        println(file, "1 bond types")
        println(file, "1 angle types\n")

        @printf file "%.8f %.8f xlo xhi\n" 0 boxSize[1,1]
        @printf file "%.8f %.8f ylo yhi\n" 0 boxSize[2,2]
        @printf file "%.8f %.8f zlo zhi\n" 0 boxSize[3,3]

        if boxSize != boxSize'
            @printf file "%.8f %.8f %.8f xy xz yz\n" boxSize[1,2] boxSize[1,3] boxSize[2,3]
        end

        println(file)

        println(file, "Masses\n")
        println(file, "1 1.00000000\n")

        println(file, "Atoms # molecular\n")

        for (atomID, (moleculeID, position)) in enumerate(atoms)
            mainImage = floorMod(position, boxSize)
            periodicIndex = floor.(boxSize\position)

            @printf file "%d %d %d" atomID moleculeID 1
            @printf file " %.8f %.8f %.8f" mainImage...
            @printf file " %d %d %d\n" periodicIndex...

            next!(p)
        end

        println(file, "\nBonds\n")

        for (bondID, bond) in enumerate(bonds)
            @printf file "%d 1 %d %d\n" bondID bond...
            next!(p)
        end

        println(file, "\nAngles\n")

        for (angleID, angle) in enumerate(angles)
            @printf file "%d 1 %d %d %d\n" angleID angle...
            next!(p)
        end

        println(file, "\n")

        finish!(p)
    end
end