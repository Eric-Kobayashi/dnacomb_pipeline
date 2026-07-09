#!/usr/bin/env bash
set -euo pipefail

DNACOMB_VERSION="${DNACOMB_VERSION:-1.0.0}"
BASE_IMAGE="${BASE_IMAGE:-mercury/dnacomb:base}"
IMAGE_TAG="${IMAGE_TAG:-dnacomb:${DNACOMB_VERSION}}"
BUILD_BASE="${BUILD_BASE:-0}"
DOCKERFILE_BASE="${DOCKERFILE_BASE:-Dockerfile-dnacomb}"
DOCKERFILE_DNACOMB="${DOCKERFILE_DNACOMB:-Dockerfile-dnacomb-addon}"
BUILD_CONTEXT="${BUILD_CONTEXT:-.}"

export DOCKER_BUILDKIT="${DOCKER_BUILDKIT:-1}"

if [[ "${BUILD_BASE}" == "1" ]]; then
  echo "Building base image from ${DOCKERFILE_BASE}: ${BASE_IMAGE}"
  docker build --pull --no-cache -f "${DOCKERFILE_BASE}" -t "${BASE_IMAGE}" "${BUILD_CONTEXT}"
else
  echo "Using existing base image: ${BASE_IMAGE}"
  docker image inspect "${BASE_IMAGE}" >/dev/null
fi

echo "Building dnacomb ${DNACOMB_VERSION} image: ${IMAGE_TAG}"
docker build --no-cache \
  -f "${DOCKERFILE_DNACOMB}" \
  --build-arg "DNACOMB_BASE_IMAGE=${BASE_IMAGE}" \
  --build-arg "DNACOMB_VERSION=${DNACOMB_VERSION}" \
  -t "${IMAGE_TAG}" \
  "${BUILD_CONTEXT}"

echo "Testing dnacomb in ${IMAGE_TAG}"
version_output="$(docker run --rm --entrypoint /bin/bash "${IMAGE_TAG}" -lc 'dnacomb --version')"
echo "${version_output}"
if [[ "${version_output}" != *" ${DNACOMB_VERSION}" ]]; then
  echo "Expected dnacomb ${DNACOMB_VERSION}, got: ${version_output}" >&2
  exit 1
fi

echo "Built ${IMAGE_TAG}"
