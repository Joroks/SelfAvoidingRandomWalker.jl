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
function saveAsLammps(self::SARW, filePath, headerTitle)
    boxSize = self.box
    chains = diff([0; self.chains])

    bonds = Int64[]
    angles = Int64[]
    
    numAtoms = 0
    for chain in chains
        for atomIndex in 1:chain
            numAtoms += 1

            atomIndex >= 2 && Base.push!(bonds, numAtoms)
            atomIndex >= 3 && Base.push!(angles, numAtoms)
        end
    end

    open(filePath, "w") do file
        desc = "Writing to LAMMPS file:"
        p = Progress(length(self.atoms)+length(bonds)+length(angles), desc)

        println(file, headerTitle)
        println(file)

        @printf(file, "%d atoms\n%d bonds\n%d angles\n\n", length(self.atoms), length(bonds), length(angles))
        println(file, "1 atom types\n1 bond types\n1 angle types\n")

        @printf(file, "%.8f %.8f %.8f avec\n", boxSize[:,1]...)
        @printf(file, "%.8f %.8f %.8f bvec\n", boxSize[:,2]...)
        @printf(file, "%.8f %.8f %.8f cvec\n", boxSize[:,3]...)

        println(file)

        println(file, "Masses\n")
        println(file, "1 1.00000000\n")

        println(file, "Atoms # molecular\n")

        moleculeID = 1
        for (atomID, pos) in enumerate(self.atoms)
            image = floor.(Int, self.box\pos)
            pos -= self.box*image

            @printf(file, "%d %d 1 %.8f %.8f %.8f %d %d %d\n", atomID, moleculeID, pos..., image...)
            moleculeID += self.chains[moleculeID] <= atomID
            next!(p)
        end

        println(file, "\nBonds\n")
        for (bondID, bond) in enumerate(bonds)
            @printf(file, "%d 1 %d %d\n", bondID, bond-1, bond)
            next!(p)
        end

        println(file, "\nAngles\n")
        for (angleID, angle) in enumerate(angles)
            @printf(file, "%d 1 %d %d %d\n", angleID, angle-2, angle-1, angle)
            next!(p)
        end

        finish!(p)
    end
end