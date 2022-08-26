# Pre-requisites

* [progressbar](https://github.com/gipert/progressbar)
    ```sh
    git clone https://github.com/gipert/progressbar
    cd progressbar

    sudo cp include/progressbar.hpp /usr/local/lib
    ```

* [cnpy](https://github.com/rogersce/cnpy)
    ```sh
    yay -S cnpy-git
    ```

* [YAM2](https://github.com/cbpark/YAM2)
    ```sh
    git clone https://github.com/cbpark/YAM2
    cd YAM2
    cd build
    cmake cmake -Dnlopt_DIR=/usr -DCMAKE_INSTALL_PREFIX=/usr/local ..
    make
    sudo make install
    ```

# Build

## Test

```sh
g++ -o bin/test test.cc -I/usr/local/include/YAM2 -L/usr/local/lib -L/usr/lib/root -lYAM2 -lnlopt -lCore 
```

## Toy

```sh
g++ -w -O3 -o bin/toy toy.cc -I/usr/include/eigen3 -L/usr/lib/root -lcnpy -lz -lYAM2 -lnlopt -lCore -std=c++17
```