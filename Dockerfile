FROM ubuntu:14.04
MAINTAINER linh lam

RUN apt-get update && apt-get install -y --force-yes --no-install-recommends \
    python \
    python-pip \
    python-dev \
    build-essential \
    zlib1g-dev \
    libcurl4-gnutls-dev \
    libssl-dev \
    pigz

WORKDIR /root
ADD https://github.com/samtools/samtools/releases/download/1.3.1/samtools-1.3.1.tar.bz2 /root
RUN tar xvf /root/samtools-1.3.1.tar.bz2
RUN make -C /root/samtools-1.3.1/htslib-1.3.1/
WORKDIR /root/samtools-1.3.1
RUN ./configure --without-curses
RUN make
RUN make install
ADD ./bamUtil /root/bamUtil
ADD ./libStatGen /root/libStatGen
ADD ./scripts /root/scripts
WORKDIR /root/libStatGen
RUN make
RUN make install
WORKDIR /root/bamUtil
RUN make
RUN make install
ADD ./run.sh /root
WORKDIR /data
ENTRYPOINT ["/root/run.sh"]
