.DEFAULT_GOAL := test

DNACOMB_VERSION ?= 1.0.0
NXF_VER ?= 26.04.6
NEXTFLOW ?= nextflow

.PHONY: test

test:
	DNACOMB_VERSION="$(DNACOMB_VERSION)" NXF_ANSI_LOG=false NXF_VER="$(NXF_VER)" "$(NEXTFLOW)" run . -with-singularity -config test/test.config
