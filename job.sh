#!/bin/bash
#SBATCH --job-name=Re5d2_PREC_ROT      # nom du job
#SBATCH --qos=qos_cpu-dev
#SBATCH -A nor@cpu
#SBATCH --ntasks=40 ##32 ###192              # Nombre total de processus MPI
####SBATCH --ntasks-per-node=40      # Nombre de processus MPI par noeud
# /!\ Attention, la ligne suivante est trompeuse mais dans le vocabulaire
# de Slurm "multithread" fait bien reference a l'hyperthreading.
#SBATCH --hint=nomultithread       # 1 processus MPI par coeur physique (pas d'hyperthreading)
#SBATCH --time=00:30:00            # Temps d execution maximum demande (HH:MM:SS)
#SBATCH --output=Re5d2_PREC_ROT%j.out  # Nom du fichier de sortie
#SBATCH --error=Re5d2_PREC_ROT%j.out   # Nom du fichier d'erreur (ici commun avec la sortie)

# on se place dans le repertoire de soumission
cd ${SLURM_SUBMIT_DIR}

# nettoyage des modules charges en interactif et herites par defaut
#module purge

# chargement des modules
#module load intel-all/19.0.4

# echo des commandes lancees
set -x
 
# exécution du code
#WORKIN=$WORK/MY_APPLICATIONS_SFEMaNS_GIT/LES_VKS_TMstraight/RUNS/Re1500_Rm150_mu50_CI1em6_3S_64F_beta_1_cdiv_1/TEST_straight_geom

date
cp $WORK/MY_APPLICATIONS_SFEMaNS_GIT/PREC_ROT_AX/EXECUTABLE/a.exe .

#cp $WORKIN/suite*I008* .
#cp $WORKIN/mesh_part_S3* .
#for FILE in suite_*I008*; do mv "$FILE" $(echo "$FILE" | sed 's/_I008//'); done
#
# execution du code
srun ./a.exe

date
