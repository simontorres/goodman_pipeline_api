FROM python:3

ENV PYTHONUNBUFFERED 0
RUN mkdir -p /app
WORKDIR /app

COPY requirements.txt /app

RUN pip install -r requirements.txt

#COPY . /api

EXPOSE 8000

#CMD ["python", "-u", "manage.py", "runserver"]