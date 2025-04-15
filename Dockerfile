# Use an official Python runtime as a parent image
# Using python:3.10-slim for a smaller image size
FROM python:3.10-slim

# Set environment variables
# Prevents Python from buffering stdout and stderr
ENV PYTHONUNBUFFERED True
# Set the working directory in the container
ENV APP_HOME /app
WORKDIR $APP_HOME

# Create a non-root user and group
# Running as non-root is a security best practice
RUN addgroup --system appgroup && adduser --system --ingroup appgroup appuser

# Install system dependencies if any were needed (none identified for this app yet)
# RUN apt-get update && apt-get install -y --no-install-recommends some-package && rm -rf /var/lib/apt/lists/*

# Copy the requirements file into the container
COPY requirements.txt .

# Install dependencies
# --no-cache-dir reduces image size
# --user installs packages for the current user (root here), but we switch later
# Alternatively, install after switching user, but might need permission adjustments
RUN pip install --no-cache-dir -r requirements.txt

# Copy the rest of the application code into the container
# Ensure you have a .dockerignore file to exclude unnecessary files (like .git, .venv, __pycache__)
COPY . .

# Change ownership of the app directory to the non-root user
RUN chown -R appuser:appgroup $APP_HOME

# Switch to the non-root user
USER appuser

# Expose the port the app runs on (standard for Cloud Run)
EXPOSE 8080

# Define the command to run the application using Gunicorn
# --workers: Adjust based on expected load and Cloud Run CPU allocation (2-4 is a common starting point)
# --bind 0.0.0.0:8080: Listen on all interfaces on the specified port
# app:app: Tells Gunicorn to run the 'app' object from the 'app.py' file
CMD ["gunicorn", "--bind", "0.0.0.0:8080", "--workers", "2", "app:app"]
