FROM soartelescope/reduction-api:1.0

ENV PYTHONUNBUFFERED 0
RUN mkdir -p /app
WORKDIR /app

COPY requirements.txt /app
RUN pip install -r requirements.txt
#RUN pip install goodman-pipeline==1.3.0.dev18

RUN pip show goodman_pipeline | grep Location | awk -v local_path=$(dirname $(which redspec)) '{ print "cp "$2"/goodman_pipeline/data/dcr-source/dcr/dcr "local_path }' | sh

COPY . /app

EXPOSE 8000

#CMD ["python", "-u", "manage.py", "runserver"]
