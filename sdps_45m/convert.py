__all__ = ["to_dems"]


# standard library
from os import PathLike
from re import search
from pathlib import Path


# dependencies
import numpy as np
from astropy.io.fits import HDUList
from dems.d2 import MS
from fmflow.fits.nro45m.functions import make_obsinfo_sam45, read_backendlog_sam45


def to_dems(
    sam45: PathLike[str] | str,
    /,
    *,
    array: str = r"A\d+",
    overwrite: bool = False,
) -> None:
    """Convert a SAM45 logging to DEMS files of each array.

    Args:
        sam45: Path of the SAM45 logging.
        array: Regular expression to select the array(s).
        overwrite: Whether to overwrite the existing DEMS files.

    """
    # read SAM45 logging using fmflow
    sam45 = Path(sam45)
    hdu_data = read_backendlog_sam45(sam45, "<")
    hdu_info = make_obsinfo_sam45(HDUList([hdu_data]))

    header = hdu_info.header
    info = hdu_info.data
    data = hdu_data.data

    # create and save DEMS of each array
    for subinfo in info:
        if not search(array, arrayid := subinfo["arrayid"]):
            continue
        else:
            subdata = data[data["arrayid"] == arrayid]

        ms = MS.new(
            # data
            data=subdata["arraydata"],
            long_name=f"{sam45.name}.{arrayid}",
            name=f"{sam45.name}.{arrayid}",
            units="dimensionless",
            # dimensions
            chan=np.arange(subinfo["chtotaln"]),
            time=np.array(subdata["starttime"], "M8[ns]") - np.timedelta64(9, "h"),
            # labels
            state=subdata["scantype"],
            scan=subdata["iline_no"],
            frequency=np.linspace(
                subinfo["rfcenter"] - subinfo["chwidth"] * (subinfo["chcenter"] - 1),
                subinfo["rfcenter"] + subinfo["chwidth"] * (subinfo["chcenter"] - 1),
                subinfo["chtotaln"],
            ),
            # telescope pointing
            lon=subdata["dant_real"][:, 0],
            lat=subdata["dant_real"][:, 1],
            lon_origin=subdata["dmap_center"][:, 2],
            lat_origin=subdata["dmap_center"][:, 5],
            # data information
            exposure=subinfo["interval"],
            interval=subinfo["interval"],
            # observation information
            object=header["object"],
            observer=header["observer"],
            telescope_name=header["telescop"],
        )

        # reassign ON/OFF states (assuming R,ON,OFF,...,ON,OFF)
        ms = ms[(ms.state == "ON") | (ms.state == "R")]
        ms.coords["state"][(ms.state == "ON") & (ms.scan.astype(int) % 2 == 0)] = "OFF"
        ms.coords["scan"][:] = np.ceil(ms.scan.astype(int) / 2).astype(int).astype(str)
        ms.coords["subscan"][:] = np.where(ms.state == "OFF", 1, 0).astype(str)

        # save DEMS as a zipped Zarr
        zarr = sam45.with_name(f"{sam45.name}.{arrayid}.zarr.zip")

        if zarr.exists() and not overwrite:
            raise FileExistsError(zarr)
        else:
            ms.to_zarr(zarr, mode="w")
