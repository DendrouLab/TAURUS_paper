#!/bin/bash 
#SBATCH -J mast_dge_trt_rmsn_int
#SBATCH -o mast_dge_trt_rmsn_int.%a.log
#SBATCH -p long
#SBATCH --cpus-per-task 6
#SBATCH -A cartography.prj.high


# If there's an error, fail the whole script
set -e -o pipefail

# Set permissions so any user in the group can 
# read/write what it's created by the script
umask 002



echo "********************************************************"
echo "* Job Details"
echo "********************************************************"
echo "Slurm job ID       : "$SLURM_ARRAY_JOB_ID
echo "Slurm task ID      : "$SLURM_ARRAY_TASK_ID
echo "Run on host      : "`hostname`
echo "Operating system : "`uname -s`
echo "Username         : "`whoami`
echo "Started at       : "`date`
echo

# Job Arguments
CELLTYPE_FILE=$1
echo "********************************************************"
echo "* Job Parameters"
echo "********************************************************"
echo "Celltypes File     : ${CELLTYPE_FILE}"
echo


# Task Arguments
CELLTYPE=$(awk -F '\t' "{if (NR==$SLURM_ARRAY_TASK_ID) print \$1}" $CELLTYPE_FILE)
INPUT_FILE=$(awk -F '\t' "{if (NR==$SLURM_ARRAY_TASK_ID) print \$2}" $CELLTYPE_FILE)
OUTPUT_PTH=$(awk -F '\t' "{if (NR==$SLURM_ARRAY_TASK_ID) print \$3}" $CELLTYPE_FILE)

echo "********************************************************"
echo "* Task Parameters"
echo "********************************************************"
echo "Celltype sce file             : ${INPUT_FILE}"
echo "Celltype                      : ${CELLTYPE}"
echo "Output Path                   : ${OUTPUT_PTH}"


# Load software modules
echo "Modules loading"
module purge
source ~/load_panpipes_public_venv

# Parameters validation
OPTIONS=""

if [[  -f "${INPUT_FILE}" ]]; then
	echo "The provided celltype input file exist: ${INPUT_FILE}."
else
    # string not provided
    echo "No celltype sce  file exists"
    exit 1
fi

if [ -d "${OUTPUT_PTH}" ]; then
    echo "Sample output directory  exists: ${OUTPUT_PTH}"
fi



echo "********************************************************"
echo "["`date`"] Running MAST DGE in R"
echo "********************************************************"
CMD1="Rscript /well/cartography/projects/analysis/302_external_oxford_datasets/TAURUS/revisions/src/001_TAURUS_DEG_MASTscript_Treatment_remisison_interaction_WITH_Site.R \
--filename ${INPUT_FILE} \
--celltype ${CELLTYPE} \
--out_dir ${OUTPUT_PTH}"

echo "Command : ${CMD1}"
echo "********************************************************"
echo
eval "${CMD1}"
echo
echo "********************************************************"
echo "["`date`"] Done"
echo "********************************************************"
exit 0
