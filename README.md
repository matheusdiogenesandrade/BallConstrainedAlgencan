# Ball Constrained Algencan 

#### Beforehand settings

Make sure you had built [scaled-algencan](https://github.com/leonardosecchin/scaled-algencan) properly.
After, you will have to build the shared library of the scaled-algencan:

```sh
cd /opt/scaled-algencan/algencan-3.1.1_scaled
make sharedlib
```

And then set the environment variables, such as below.

```sh
export SCALED_ALGENCAN_LIB=/opt/scaled-algencan/algencan-3.1.1_scaled/lib/
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$SCALED_ALGENCAN_LIB
```

## _CMake_

When using CMake, there is no need for setting neither `SCALED_ALGENCAN_LIB` nor `LD_LIBRARY_PATH`,
since the `FindSCALED_ALGENCAN` will load it automatically.

```sh
cd cmake/
mkdir build
cd build
cmake ..
make
./main
```


