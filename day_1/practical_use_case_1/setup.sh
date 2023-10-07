#!/bin/bash


# Linux Intel (x86_64):
curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xvj bin/micromamba

# Linux/bash:
./bin/micromamba shell init -s bash -p ~/micromamba  # this writes to your .bashrc file
# sourcing the bashrc file incorporates the changes into the running session.
# better yet, restart your terminal!
source ~/.bashrc

# activates the base environment
microbamba activate

# defines channels
micromamba config append channels bioconda
micromamba config append channels conda-forge
micromamba config append channels nodefaults
micromamba config set channel_priority strict

microbamba install --yes bcftools