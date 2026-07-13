DNACOMB_VERSION ?= 1.0.0
NXF_VER ?= 26.04.6
NEXTFLOW ?= nextflow

.PHONY: init release smoketest

init:
	./init-dev.sh

release: init
	./release-dnacomb-images.sh release

smoketest:
	DNACOMB_VERSION="$(DNACOMB_VERSION)" NXF_ANSI_LOG=false NXF_VER="$(NXF_VER)" "$(NEXTFLOW)" run .. -with-singularity -config data/smoke.config -work-dir data/work
