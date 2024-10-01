#!/bin/bash

# Define paths to your data and scripts
MAIN_DIR="$(pwd)"
DATA_DIR="${MAIN_DIR}/data/fastqs"
QC_TRIM_SCRIPT="qc_and_trim.sh"
ALIGN_SCRIPT="align.sh"
DOCKER_IMAGE="rnaseq"

# Check if Docker image exists; if not, build it
if [[ "$(docker images -q ${DOCKER_IMAGE} 2> /dev/null)" == "" ]]; then
  echo "Docker image ${DOCKER_IMAGE} not found. Building the image..."
  docker build -t ${DOCKER_IMAGE} ..  || { echo "Docker build failed. Exiting..."; exit 1; }
fi

# Step 1: Run FastQC and trimming
echo "Running FastQC and trimming..."
docker run --rm -v "${MAIN_DIR}:/scripts" \
               -v "${DATA_DIR}:/data" \
               --user $(id -u):$(id -g) \
               ${DOCKER_IMAGE} /bin/bash scripts/${QC_TRIM_SCRIPT} size_project data/raw/ data/trimmed data/fastqc

# if [ $? -ne 0 ]; then
#   echo "FastQC and trimming failed. Check the log output above for more details. Exiting..."
#   exit 1
# fi

# # Step 2: Run STAR alignment
# echo "Running STAR alignment..."
# docker run --rm -v "${ALIGN_SCRIPT}:/align.sh" \
#                -v "${DATA_DIR}:/data" \
#                ${DOCKER_IMAGE} bash /align.sh 

# if [ $? -ne 0 ]; then
#   echo "STAR alignment failed. Check the log output above for more details. Exiting..."
#   exit 1
# fi

echo "RNA-seq pipeline completed successfully."
