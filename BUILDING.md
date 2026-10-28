# Building and checking winratioIC on Windows

The package source tree is the `winratioIC` folder.

## 1. Open a terminal in the folder that contains `winratioIC`

For example:

```bat
cd "C:\path\to\simulation\software"
```

## 2. Build

```bat
R CMD build winratioIC
```

This should create:

```text
winratioIC_1.0.0.tar.gz
```

## 3. Check

```bat
R CMD check winratioIC_1.0.0.tar.gz
```

The desired final line is `Status: OK`.

## 4. Install the built source package

In a fresh R session:

```r
install.packages("winratioIC_1.0.0.tar.gz", repos = NULL, type = "source")
library(winratioIC)
?winratio_ic
```

Then run the examples in `README.md`.
