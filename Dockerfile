FROM bentoml/model-server:0.11.0-py37
MAINTAINER ersilia

#RUN python -m pip install --upgrade pip==20.2.3

RUN pip install rdkit-pypi==2021.3.5.1
#RUN conda install -c conda-forge rdkit::rdkit==2020.03.1
RUN pip install numpy==1.19.2
RUN pip install hypopt==1.0.9
RUN pip install argparse==1.4.0
RUN pip install tqdm==4.49.0
RUN pip install flask
RUN pip install cairosvg==2.4.2
RUN pip install requests
RUN pip install pubchempy==1.0.4
RUN pip install func_timeout==4.3.5
RUN pip install xgboost==1.0.2
RUN pip install scikit-learn==0.22.1
RUN pip install pandas==1.1.2


WORKDIR /repo
COPY . /repo
