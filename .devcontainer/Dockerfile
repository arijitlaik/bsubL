# Dockerfile
FROM julia:latest

# Install necessary packages
RUN apt-get update && apt-get install -y git sudo zsh curl

# Create a non-root user to use if preferred - see https://aka.ms/vscode-remote/containers/non-root-user.
RUN useradd -m codespace && \
    echo codespace:codespace | chpasswd && \
    echo "codespace ALL=(root) NOPASSWD:ALL" > /etc/sudoers.d/codespace && \
    chmod 0440 /etc/sudoers.d/codespace && \
    usermod --shell /bin/bash codespace

# Install Oh My Zsh for a better terminal experience
RUN curl -Lo install.sh https://raw.githubusercontent.com/ohmyzsh/ohmyzsh/master/tools/install.sh && \
    chmod +x install.sh && \
    sh install.sh --unattended && \
    rm install.sh

# Install Julia packages
RUN julia -e 'using Pkg; Pkg.add(["Plots", "PolygonOps"])'

# Switch back to dialog for any ad-hoc use of apt-get
ENV DEBIAN_FRONTEND=dialog

# Set the default shell to bash rather than sh
ENV SHELL /bin/bash

USER codespace
