# Dockerfile
FROM python:3.9

# Install required Python packages
RUN pip install --no-cache-dir matplotlib pandas seaborn

# Set the working directory
WORKDIR /usr/src/app

# Copy the Python script to the container
COPY alignment_stats_plot.py .

# Command to run the Python script
CMD ["python", "alignment_stats_plot.py"]

