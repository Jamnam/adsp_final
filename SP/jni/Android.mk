LOCAL_PATH := $(call my-dir)

include $(CLEAR_VARS)

LOCAL_MODULE    := libSpeechProcessing
LOCAL_SRC_FILES := RandomForestClassifier.c  SubBandFeatures.c vad.c Timer.c SpeechProcessing.c periodogram.c  Traindata.c logMMSE.c 
LOCAL_CFLAGS := -O3

LOCAL_LDLIBS := -llog -landroid

include $(BUILD_SHARED_LIBRARY)
