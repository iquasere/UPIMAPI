FROM continuumio/miniconda3

RUN git clone https://github.com/iquasere/UPIMAPI.git \
&& conda env update --file UPIMAPI/cicd/environment.yml --name base \
&& bash UPIMAPI/cicd/ci_build.sh \
&& conda clean --all -y

CMD [ "python", "bin/upimapi.py" ]