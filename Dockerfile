FROM bentoml/model-server:0.11.0-py38
MAINTAINER ersilia

RUN pip install rdkit==2023.9.5
RUN pip install scikit-learn==0.22.1
RUN pip install numpy==1.19.5
RUN pip install pandas==1.1.2
RUN pip install hypopt==1.0.9
RUN pip install xgboost==1.0.2

WORKDIR /repo
COPY . /repo
