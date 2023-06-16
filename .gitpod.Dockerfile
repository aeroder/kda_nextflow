FROM gitpod/workspace-full
#FROM rocker/r-ubuntu:18.04

RUN brew install R
RUN Rscript -e "install.packages(c('class','cluster','rpart','lattice','box'), repos = 'http://cran.us.r-project.org')"