#Used for common JASP qmake build settings

#Jasp-R-Interface
JASP_R_INTERFACE_TARGET = JASP-R-Interface

JASP_R_INTERFACE_MAJOR_VERSION =  8  # Interface changes
JASP_R_INTERFACE_MINOR_VERSION =  1 # Code changes

JASP_R_INTERFACE_NAME = $$JASP_R_INTERFACE_TARGET$$JASP_R_INTERFACE_MAJOR_VERSION'.'$$JASP_R_INTERFACE_MINOR_VERSION

#R settings
CURRENT_R_VERSION = 3.6
DEFINES += "CURRENT_R_VERSION=\"$$CURRENT_R_VERSION\""

#JASP Version
JASP_VERSION_MAJOR      = 0
JASP_VERSION_MINOR      = 11
JASP_VERSION_REVISION   = 1
JASP_VERSION_BUILD      = 0 #Should be ignored because the code handling it is buggy as hell (aka https://www.youtube.com/watch?v=otCpCn0l4Wo )

DEFINES +=    "JASP_VERSION_MAJOR=$$JASP_VERSION_MAJOR"
DEFINES +=    "JASP_VERSION_MINOR=$$JASP_VERSION_MINOR"
DEFINES +=    "JASP_VERSION_BUILD=$$JASP_VERSION_BUILD"
DEFINES +=    "JASP_VERSION_REVISION=$$JASP_VERSION_REVISION"

equals(ROOT_DIR, "") {
    ROOT_DIR =  $$OUT_PWD/../..
}

JASP_REQUIRED_FILES = $$ROOT_DIR/jasp-required-files

message(using JASP_REQUIRED_FILES of $$JASP_REQUIRED_FILES)

INCLUDEPATH += $$JASP_REQUIRED_FILES/boost_1_71_0/

BUILDING_JASP_ENGINE=false
DEFINES += PRINT_ENGINE_MESSAGES

#macx | windows | exists(/app/lib/*) {
#	message(using libjson static)
  DEFINES += JASP_LIBJSON_STATIC # lets just always use libjson-static, they keep moving the include files...
#} else {
#    linux {
#        message(using libjson from distro and pkgconfig)
#        QT_CONFIG -= no-pkg-config
#        CONFIG += link_pkgconfig
#        PKGCONFIG += jsoncpp
#        LIBS += -ljsoncpp
#
#        CONFIG(debug, debug|release) {  DEFINES+=JASP_DEBUG }
#    }
#}

exists(/app/lib/*) {
  linux:  DEFINES += FLATPAK_USED
} else {
  linux:	CONFIG(debug, debug|release)  {
    DEFINES += JASP_DEBUG
    DEFINES += LINUX_NOT_FLATPAK
  }
}
macx | windows { CONFIG(debug, debug|release) {  DEFINES += JASP_DEBUG } }

windows {
	message(QT_ARCH $$QT_ARCH)
	contains(QT_ARCH, i386) {
		ARCH = i386
		BOOST_ARCH = x32
	} else {
		ARCH = x64
		BOOST_ARCH = x64
	}

	CONFIG(ReleaseBuild)    { BOOST_POSTFIX = -vc141-mt-$${BOOST_ARCH}-1_71 }
	CONFIG(DebugBuild)      { BOOST_POSTFIX = -vc141-mt-gd-$${BOOST_ARCH}-1_71 }

	QMAKE_CXXFLAGS	+= -DBOOST_USE_WINDOWS_H -DNOMINMAX -DBOOST_INTERPROCESS_BOOTSTAMP_IS_SESSION_MANAGER_BASED
}

unix: QMAKE_CXXFLAGS += -Werror=return-type

#want to use JASPTIMER_* ? set JASPTIMER_USED to true, run qmake and rebuild the objects that use these macros (or just rebuild everything to be sure)
JASPTIMER_USED = false

$$JASPTIMER_USED {
    DEFINES += PROFILE_JASP
}

exists(/app/lib/*)	{
  INSTALLPATH = /app/bin
 } else	{
  INSTALLPATH = /usr/bin
}

DEFINES += QT_NO_FOREACH #Come on Qt we can just use the nice new ranged for from c++11 and higher, we dont need your help!

macx {
	QMAKE_CXXFLAGS_WARN_ON	+= -Wno-unused-parameter -Wno-unused-local-typedef
	QMAKE_CXXFLAGS			+= -Wno-c++11-extensions -Wno-c++11-long-long -Wno-c++11-extra-semi -stdlib=libc++
}
