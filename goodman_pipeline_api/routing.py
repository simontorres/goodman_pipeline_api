from django.conf.urls import url
from channels.routing import ProtocolTypeRouter, URLRouter
from channels.auth import AuthMiddlewareStack

from api.consumer import ApiMessageConsumer

live_reduction = ProtocolTypeRouter({
    "websocket": AuthMiddlewareStack(
        URLRouter([
            url('ws/messsages', ApiMessageConsumer),
        ])
    )
})
