#  qmake project file for the gemcalc example (part of GEMS3K standalone code)
# (c) 2012-2019 GEMS Developer Team
 
TEMPLATE = app
LANGUAGE = C++
TARGET = gemcalc
VERSION = 3.4.6

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

SOURCES      +=   main.cpp

include($$GEMS3K_CPP/gems3k.pri) 
