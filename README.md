# sdps-45m
Tools of sdps for the Nobeyama 45m Telescope

## Installation

```shell
pip install sdps-45m
```

## Usage

### Convert SAM45 logging to DEMS files

```shell
sdps-45m convert /path/to/sam45 [--overwrite]
```

This will create `/path/to/sam45.[A1-A16].zarr.zip`.
