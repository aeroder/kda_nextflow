FROM rocker/r-ubuntu:18.04

RUN Rscript -e "install.packages(c('class','cluster','rpart','lattice','box','glue'), repos = 'http://cran.us.r-project.org')"