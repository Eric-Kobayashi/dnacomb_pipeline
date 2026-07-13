DNACOMB_VERSION ?= 1.0.0
NXF_VER ?= 26.04.6
NEXTFLOW ?= nextflow

.PHONY: init release

init:
	./init-dev.sh

release: init
	./release-dnacomb-images.sh release
