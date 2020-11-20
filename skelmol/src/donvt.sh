

## do nvt md preparation from pdb
## type 15 for the force field

name="1amy0"

gmx pdb2gmx -f ${name}.pdb -o ${name}.gro -water spce -p ${name}.top
gmx editconf -f ${name}.gro -o ${name}_newbox.gro -c -d 1.0 -bt cubic
gmx solvate -cp ${name}_newbox.gro -cs spc216.gro -o ${name}_solv.gro -p ${name}.top
gmx grompp -f ions.mdp -c ${name}_solv.gro -p ${name}.top -o ions.tpr
gmx genion -s ions.tpr -o ${name}_solv_ions.gro -p ${name}.top -pname NA -nname CL -neutral
gmx grompp -f minim.mdp -c ${name}_solv_ions.gro -p ${name}.top -o ${name}_em.tpr
gmx mdrun -v -deffnm ${name}_em
gmx grompp -f nvt300.mdp -c ${name}_em.gro -p ${name}.top -o ${name}_nvt300.tpr
gmx mdrun -v -deffnm ${name}_nvt300


##

gmx grompp -f npt300.mdp -c ${name}_nvt300.gro -t ${name}_nvt300.cpt -p ${name}.top -o ${name}_npt300.tpr
gmx mdrun -v -deffnm ${name}_npt300

