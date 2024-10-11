# Use an official Python base image with a suitable version for your needs
FROM debian:bullseye-slim

# Set the working directory inside the container
WORKDIR /

# Install system dependencies including build tools
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    wget \
    curl \
    unzip \
    tar \
    build-essential \
    zlib1g-dev \
    gcc \
    g++ \
    libz-dev \
    openjdk-11-jre-headless \
    make \
    python3 \
    python3-pip \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# Install FastQC
RUN wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip \
    && unzip fastqc_v0.11.9.zip \
    && rm fastqc_v0.11.9.zip \
    && chmod +x FastQC/fastqc \
    && mv FastQC /usr/local/bin/

# Install MultiQC
RUN pip install multiqc==1.9.0

# Install Trim Galore
RUN pip install cutadapt
RUN wget https://github.com/FelixKrueger/TrimGalore/archive/0.6.7.tar.gz \
    && tar -xzvf 0.6.7.tar.gz \
    && rm 0.6.7.tar.gz \
    && cp TrimGalore-0.6.7/trim_galore /usr/local/bin/ \
    && chmod +x /usr/local/bin/trim_galore \
    && rm -rf TrimGalore-0.6.7

# Download and install STAR
RUN wget https://github.com/alexdobin/STAR/archive/refs/tags/2.7.10a.tar.gz \
    && tar -xzvf 2.7.10a.tar.gz \
    && rm 2.7.10a.tar.gz \
    && cd STAR-2.7.10a/bin/Linux_x86_64_static \
    && make STAR \
    && cp STAR /usr/local/bin/ \
    && rm -rf ../../../STAR-2.7.10a

# Copy the requirements.txt file into the container
COPY requirements.txt /usr/local/bin

# Install Python packages from the requirements.txt file
RUN pip install -r /usr/local/bin/requirements.txt

# Set softwares binary location to PATH
ENV PATH="/usr/local/bin/FastQC:${PATH}"
ENV PATH="/usr/local/bin/TrimGalore:${PATH}"
ENV PATH="/usr/local/bin/STAR:${PATH}"

# Set the JAVA_OPTS environment variable to run in headless mode
ENV JAVA_OPTS="-Djava.awt.headless=true"

# Set default command for container
CMD ["bash"]
