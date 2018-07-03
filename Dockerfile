# copyright 2017-2018 Regents of the University of California and the Broad Institute. All rights reserved.

FROM r-base:3.5.0

RUN mkdir /build

RUN apt-get update && apt-get upgrade --yes && \
    apt-get install -t unstable libmariadbclient-dev  --yes && \
    apt-get install -t unstable libssl-dev  --yes && \
    apt-get install libxml2-dev --yes && \
    apt-get install libcurl4-gnutls-dev --yes && \
    apt-get install mesa-common-dev --yes && \ 
    apt-get update && apt-get install -y --no-install-recommends apt-utils && \
    apt-get install aptitude -y && \
    apt-get install libxml2-dev -y && \
    aptitude install libglib2.0-dev -y && \
    apt-get install libcairo2-dev -y && \
    apt-get install  xvfb xauth xfonts-base libxt-dev -y && \
    rm -rf /var/lib/apt/lists/*

COPY sources.list /etc/apt/sources.list
COPY Rprofile.gp.site ~/.Rprofile
COPY Rprofile.gp.site /usr/lib/R/etc/Rprofile.site
COPY install_stuff.R /build/source/install_stuff.R
RUN Rscript /build/source/install_stuff.R

ENV PATH "$PATH:/usr/local/bin/sctk"
COPY src/* /usr/local/bin/sctk/

 
CMD ["Rscript", "/usr/local/bin/cogaps/run_gp_tutorial_module.R" ]

