#  qmake project file for the gems3_server example (part of GEMS3K standalone code)
# (c) 2020 GEMS Developer Team
 
TEMPLATE = app
LANGUAGE = C++
TARGET = gems3_client
VERSION = 0.1

CONFIG -= qt
CONFIG -= warn_on
CONFIG += debug
#CONFIG += windows
CONFIG += console


DEFINES += IPMGEMPLUGIN
#DEFINES += NODEARRAYLEVEL
DEFINES += NOPARTICLEARRAY

!win32 {
  DEFINES += __unix
}

GEMS3K_CPP = ../GEMS3K
GEMS3K_H   = $$GEMS3K_CPP

DEPENDPATH +=
DEPENDPATH += .
DEPENDPATH += $$GEMS3K_H

INCLUDEPATH += 
INCLUDEPATH += .
INCLUDEPATH += $$GEMS3K_H

QMAKE_LFLAGS +=
QMAKE_CXXFLAGS += -Wall -Wno-unused
OBJECTS_DIR = obj

SOURCES      +=   client_main.cpp

LIBS += -lzmq

include($$GEMS3K_CPP/gems3k.pri) 
