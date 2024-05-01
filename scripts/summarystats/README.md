# Producing a standard set of image summary statistics

This script produces a standard set of image summary statistics for a list of images. To run the script first call `setup.jl` from the command line

```bash
julia setup.jl
```

After the setup is complete you can run the main script `main.jl` to produce the set of summary statistics. To run the script call

```bash
julia -p NCORES main.jl filelist output.csv
```
where `NCORES` is the number of cores, `filelist` is a text file containing paths, the list of images to process, and `output.csv` is the name of the output file. In addition, there are several other options and flags

```
# Options

- `-c, --code=<string>`: The file containing the code used to generate the images. If nothing will ignore.
- `-f, --fevals=<int>`: The number of evaluations allowed when extracting params
- `-s, --stride=<int>`: The checkpointing stride
- `-o, --order=<int>`: The order of the ring parameters
- `-b, --blur=<float>`: Gaussian FWHM blurring kernel to apply to the images in μas

# Flags

- `-r, --regrid`: Regrid the images before extracting
- `--restart`: Restart the extraction process from before
```
