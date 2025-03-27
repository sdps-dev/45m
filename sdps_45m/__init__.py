__all__ = ["convert"]


# dependencies
from fire import Fire
from . import convert


# command line interface
def main() -> None:
    Fire({"convert": convert.to_dems})


if __name__ == "__main__":
    main()
