# sdps-45m

[![Release](https://img.shields.io/pypi/v/sdps-45m?label=Release&color=cornflowerblue&style=flat-square)](https://pypi.org/project/sdps-45m/)
[![Python](https://img.shields.io/pypi/pyversions/sdps-45m?label=Python&color=cornflowerblue&style=flat-square)](https://pypi.org/project/sdps-45m/)
[![Downloads](https://img.shields.io/pypi/dm/sdps-45m?label=Downloads&color=cornflowerblue&style=flat-square)](https://pepy.tech/project/sdps-45m)
[![Tests](https://img.shields.io/github/actions/workflow/status/sdps-dev/45m/tests.yaml?label=Tests&style=flat-square)](https://github.com/sdps-dev/45m/actions)

Tools of sdps for the Nobeyama 45m Telescope

## Installation

```shell
pip install sdps-45m
```

## Usage

### Convert SAM45 log to DEMS files

```shell
sdps-45m convert /path/to/sam45 [--overwrite]
```

This will create `/path/to/sam45.[A1-A16].zarr.zip`.
