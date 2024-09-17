# Dockerfile
FROM python:3.9

# Install required Python packages
RUN pip install --no-cache-dir matplotlib pandas seaborn

# Install bwa and samtools
RUN apt-get update && \
    apt-get install -y bwa samtools

# Set the working directory
WORKDIR /usr/src/app

# Copy the Python script to the container
COPY plot_alignment.py .
COPY calculate_percentage.py .

# Command to run the Python script
CMD ["python", "plot_alignment.py"]

# Command to run the Python script
CMD ["python", "calculate_percentage.py"]

