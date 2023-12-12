#!/bin/bash
# Check if this directory contains: PS2_1.mol2/.pdb/.itp/.prm/.top  binary_1.pdb ps2_1_ini.pdb
# Check if both small molecule name and $i value are provided as arguments
set -e

if [ -z "$1" ] || [ -z "$2" ]; then
    echo 'Please provide both the small molecule name and the value for $i as arguments.'
    exit 1
fi

# Assign the small molecule name and $i value to variables
molecule_name="$1"
i="$2"

if [ ! -d "charmm36-jul2021.ff" ]; then
        ln -s ~/charmm36-jul2021.ff ./
fi
gmx editconf -f "${molecule_name}_${i}_ini.pdb" -o "${molecule_name}_${i}.gro"
sed -i "s/ZN    ZN/ZN   ZN2/g" binary_$i.pdb
cp /data/apps/sbgrid/x86_64-linux/gromacs/2022.1_cu11.5.2/share/gromacs/top/residuetypes.dat ./
sed -i "s/ZN2/ZN/g" charmm36-jul2021.ff/ions.itp
sed -i '$aZN2     Ion\
SOD     Ion\
CLA     Ion' residuetypes.dat
gmx pdb2gmx -f binary_$i.pdb -o binary_$i.gro -ignh -ff charmm36-jul2021 -water tip3p
sleep 5
cp binary_$i.gro ternary_$i.gro
tail -n +3 "${molecule_name}_$i.gro" | head -n -1 > lig1coord
sed -i "/^$(tail -n2 ternary_$i.gro | head -n1)/ r lig1coord" ternary_$i.gro
totalatom=$(expr $(sed -n 2p ternary_$i.gro) + $(sed -n 2p "${molecule_name}_$i.gro"))
sed -i "2c $totalatom" ternary_$i.gro
sed -i 's/; Include chain topologies/; Include ligand parameters\
#include "'"$molecule_name"'_'"$i"'.prm"\
\
; Include chain topologies/g' topol.top
sed -i 's/; Include water topology/; Include ligand topology\
#include "'"$molecule_name"'_'"$i"'.itp"\
\
; Include water topology/g' topol.top
sed -i '$a'"$molecule_name"'             1' topol.top
cp -n ~/mdp_protein_ligand/*.mdp ./
sed -i "s/mz1/$molecule_name/g" *.mdp
gmx editconf -f ternary_$i.gro -o newbox.gro -bt cubic -d 1.0
gmx solvate -cp newbox.gro -cs spc216.gro -p topol.top -o solv.gro
gmx editconf -f solv.gro -o solv.pdb
gmx grompp -f ions.mdp -c solv.gro -p topol.top -o ions.tpr
echo "SOL" | gmx genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -neutral
gmx grompp -f em.mdp -c solv_ions.gro -p topol.top -o em.tpr
sleep 5
gmx_mpi mdrun -v -s em.tpr -deffnm em -pin on
echo non-Water | gmx trjconv -f em.gro -s em.tpr -o em.pdb
( echo '0&!aH*'; echo q ) | gmx make_ndx -f "${molecule_name}_$i.gro" -o index_"$molecule_name"_$i.ndx
echo 3 | gmx genrestr -f "${molecule_name}_$i.gro" -n index_"$molecule_name"_$i.ndx -o posre_"$molecule_name".itp -fc 1000 1000 1000
sed -i 's/; Include water topology/; Ligand position restraints\
#ifdef POSRES\
#include "posre_'"$molecule_name"'.itp"\
#endif\
\
; Include water topology/g' topol.top
echo q | gmx make_ndx -f em.gro -o index.ndx
gmx select -f em.gro -s em.tpr -on index_protein_"$molecule_name".ndx -select 'resname '"$molecule_name"' or group 1'
sed -i 's/resname_'"$molecule_name"'_or_group_1/Protein_'"$molecule_name"'/g' index_protein_"$molecule_name".ndx
cat index_protein_"$molecule_name".ndx >> index.ndx
gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -n index.ndx -o nvt.tpr
sleep 5
export OMP_NUM_THREADS=128
gmx_mpi mdrun -v -s nvt.tpr -deffnm nvt -pin on
sed -i 's/Berendsen/C-rescale/g' npt.mdp
gmx grompp -f npt.mdp -c nvt.gro -t nvt.cpt -r nvt.gro -p topol.top -n index.ndx -o npt.tpr
sleep 5
gmx_mpi mdrun -v -s npt.tpr -deffnm npt
sed -i -e 's/500000 /50000000 /g' -e 's/1000 ps/100,000 ps/g' -e 's/1 ns/100 ns/g' md.mdp
sed -i -e 's/5000 /50000 /g' -e 's/10.0 ps/100 ps/g' md.mdp
gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -n index.ndx -o md_"$molecule_name".tpr
