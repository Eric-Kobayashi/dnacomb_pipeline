#!/usr/bin/env bash
set -euo pipefail

CRATE_NAME="${CRATE_NAME:-dnacomb}"
DOCKERHUB_USER="${DOCKERHUB_USER:-$(docker info 2>/dev/null | awk -F': ' '/^ Username:/ {print $2; exit}')}"
IMAGE_NAME="${IMAGE_NAME:-dnacomb}"
IMAGE_REPO="${IMAGE_REPO:-${DOCKERHUB_USER}/${IMAGE_NAME}}"
BASE_IMAGE="mercury/dnacomb:base"
TOKEN_FILE="${TOKEN_FILE:-docker.token}"
CRATES_USER_AGENT="${CRATES_USER_AGENT:-dnacomb-image-release/1.0}"

if [[ -z "${DOCKERHUB_USER}" ]]; then
  echo "Set DOCKERHUB_USER or run docker login first." >&2
  exit 1
fi

if [[ ! -s "${TOKEN_FILE}" ]]; then
  echo "Missing DockerHub token file: ${TOKEN_FILE}" >&2
  exit 1
fi

tmp_docker_config=""
if [[ -z "${DOCKER_CONFIG:-}" ]]; then
  tmp_docker_config="$(mktemp -d)"
  export DOCKER_CONFIG="${tmp_docker_config}"
  trap 'rm -rf "${tmp_docker_config}"' EXIT
fi

docker login --username "${DOCKERHUB_USER}" --password-stdin < "${TOKEN_FILE}" >/dev/null

remote_tag_exists() {
  docker manifest inspect "docker.io/${IMAGE_REPO}:$1" >/dev/null 2>&1
}

crate_versions() {
  curl -fsSL -H "User-Agent: ${CRATES_USER_AGENT}" "https://crates.io/api/v1/crates/${CRATE_NAME}" |
    jq -r '.versions[] | select(.yanked == false) | .num' |
    sort -Vr
}

release_versions() {
  echo "Base image: ${BASE_IMAGE}"
  echo "Release repo: ${IMAGE_REPO}"

  versions="$(crate_versions)"
  if [[ -z "${versions}" ]]; then
    echo "No ${CRATE_NAME} versions found." >&2
    exit 1
  fi

  while IFS= read -r version; do
    if remote_tag_exists "${version}"; then
      echo "Exists: ${IMAGE_REPO}:${version}"
      continue
    fi

    echo "Missing: ${IMAGE_REPO}:${version}"
    DNACOMB_VERSION="${version}" \
      BASE_IMAGE="${BASE_IMAGE}" \
      IMAGE_TAG="${IMAGE_REPO}:${version}" \
      ./build-dnacomb-image.sh
    docker push "${IMAGE_REPO}:${version}"
  done <<< "${versions}"
}

case "${1:-release}" in
  release)
    release_versions
    ;;
  *)
    echo "Usage: $0 [release]" >&2
    exit 2
    ;;
esac
