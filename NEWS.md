# winratioIC 1.0.0

* Companion-software release for the revised methodology.
* Supports interval/interval and right/interval prioritized outcomes.
* Adds user-selectable lower integration bound `tau_l`, with `0 <= tau_l < tau`.
* For interval/interval, both prioritized components use `[tau_l, tau]`.
* For right/interval, the primary component remains on `[0, tau]` and the secondary interval-censored component uses `[tau_l, tau]`.
* `tau_l = 0` computes the untrimmed exploratory target.
* Includes the revised finite-horizon dual-interval win functional and revised mixed-censoring quadrature.
* Includes weighted-bootstrap inference and fit diagnostics.
