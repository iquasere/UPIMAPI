FROM continuumio/miniconda3

RUN git clone https://github.com/iquasere/UPIMAPI.git \
&& conda install -c conda-forge mamba libarchive=3.6.2=h039dbb9_1 \
&& mamba env update --file UPIMAPI/cicd/environment.yml --name base \
&& bash UPIMAPI/cicd/ci_build.sh \
&& conda clean --all -y

CMD [ "python", "bin/upimapi.py" ]