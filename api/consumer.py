from channels.generic.websocket import WebsocketConsumer


class ApiMessageConsumer(WebsocketConsumer):

    def connect(self):
        self.accept()

    def receive(self, text_data=None, bytes_data=None):
        self.send("hello from api")

    def disconnect(self, code):
        pass
