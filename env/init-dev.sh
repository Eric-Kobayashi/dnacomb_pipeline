#!/usr/bin/env bash
set -euo pipefail

DOCKERHUB_USER="${DOCKERHUB_USER:-$(docker info 2>/dev/null | awk -F': ' '/^ Username:/ {print $2; exit}')}"
TOKEN_FILE="${TOKEN_FILE:-docker.token}"
PIXI_BIN="${PIXI_BIN:-}"

need_cmd() {
  if ! command -v "$1" >/dev/null 2>&1; then
    echo "Missing required command: $1" >&2
    exit 1
  fi
}

install_pixi() {
  need_cmd curl
  curl -fsSL https://pixi.sh/install.sh | sh
}

find_pixi() {
  if [[ -n "${PIXI_BIN}" && -x "${PIXI_BIN}" ]]; then
    return 0
  fi
  if command -v pixi >/dev/null 2>&1; then
    PIXI_BIN="$(command -v pixi)"
    return 0
  fi
  if [[ -x "${HOME}/.pixi/bin/pixi" ]]; then
    PIXI_BIN="${HOME}/.pixi/bin/pixi"
    return 0
  fi
  return 1
}

ensure_pixi() {
  if ! find_pixi; then
    echo "Pixi not found; installing."
    install_pixi
    find_pixi
  fi

  echo "Using Pixi: ${PIXI_BIN}"
  "${PIXI_BIN}" install
}

ensure_docker_login() {
  if docker info 2>/dev/null | awk -F': ' '/^ Username:/ {found=1} END {exit !found}'; then
    echo "DockerHub login present."
    return 0
  fi

  echo "DockerHub login missing."
  if [[ -s "${TOKEN_FILE}" ]]; then
    read -r -p "DockerHub username${DOCKERHUB_USER:+ [${DOCKERHUB_USER}]}: " entered_user
    DOCKERHUB_USER="${entered_user:-${DOCKERHUB_USER}}"
    if [[ -z "${DOCKERHUB_USER}" ]]; then
      echo "DockerHub username required." >&2
      exit 1
    fi
    docker login --username "${DOCKERHUB_USER}" --password-stdin < "${TOKEN_FILE}"
  else
    echo "No ${TOKEN_FILE}; starting interactive docker login."
    docker login
  fi
}

need_cmd bash
need_cmd curl
need_cmd docker
need_cmd jq
need_cmd make

docker info >/dev/null
ensure_pixi
ensure_docker_login

echo "Init OK."
