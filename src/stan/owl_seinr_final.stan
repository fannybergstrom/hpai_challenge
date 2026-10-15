// ============================================================
// Spatial SEINR transmission model for HPAI on poultry premises
//
// Observed per premises: t_suspicion, t_notification, t_cull_start,
//                        t_cull_complete, t_prev_cull, prev_cull_state
// Latent per premises:   infection day (S -> E), infectious day (E -> I)
//
// Likelihood cases:
//   notified                 infection and infectious days marginalised over
//                            their feasible windows
//   preventively culled in I infectious period right-censored at the cull
//   preventively culled in E latent period right-censored at the cull
//   otherwise                survival (escape from infection)
//
// active[i, t] > 0 is the relative flock size (0 = empty or culled).
//
// Duration distributions (discrete):
//   E: NegBin(r = 4, p = 0.5), mean 4 days (sigma_r, sigma_p)
//   I: geometric(1/mu)
//   N (notification to cull start): geometric(1/rho_n)
//   culling: lognormal
// ============================================================

functions {
  // Log probability that a geometric(p) r.v. equals k (k >= 1)
  real log_geom(int k, real p) {
    return log1m(p) * (k - 1) + log(p);
  }

  // rho_lo before t_change, rho_hi from t_change onward
  real rho_n_step(int t, int t_change, real rho_lo, real rho_hi) {
    return (t < t_change) ? rho_lo : rho_hi;
  }

  // Log PMF of NegBin(r, p) at k >= 0 (failures before r successes)
  real log_negbin_duration(int k, real r, real p) {
    return lgamma(r + k) - lgamma(r) - lgamma(k + 1)
           + r * log(p) + k * log1m(p);
  }
}

data {
  // ---- Dimensions ----
  int<lower=1> F;   // total farms
  int<lower=1> T;   // fitting window length (days)
  int<lower=1> K;   // number production types

  // ---- Farm attributes ----
  array[F] int<lower=1, upper=K> prod_type;
  array[F] int<lower=0, upper=1> in_hrz;  // indicator for high risk zone

  // ---- Premises activity (relative flock size; 0 = no birds) ----
  matrix<lower=0>[F, T] active;

  // ---- Sparse spatial neighbour list (CSR format, pairs within 50 km) ----
  // Kernel K(d) = exp(-delta * d)
  int<lower=0>           n_neighbors;
  array[n_neighbors] int nbr_i;        // receptor farm index (1-based)
  vector[n_neighbors]    nbr_dist_km;  // distance in km
  array[F+1] int         j_ptr;        // CSR row pointer (1-based, length F+1)

  // ---- Observations (0 = not observed / farm never infected) ----
  // Days are indexed from the start of the fitting period
  array[F] int<lower=0> t_suspicion;       // day of suspicion report
  array[F] int<lower=0> t_notification;    // day confirmed (I -> N)
  array[F] int<lower=0> t_cull_start;      // day culling began (N -> R start)
  array[F] int<lower=0> t_cull_complete;   // day culling finished

  // ---- Preventive cull ----
  array[F] int<lower=0> t_prev_cull;        // day preventively culled (0 if none)
  // Compartment at preventive cull: 0 = S, 2 = E, 3 = I
  array[F] int<lower=0, upper=3> prev_cull_state;

  // ---- Confinement ----
  array[K] int confined_from;    // day confinement started per type (T+1 if never)
  array[K] int confined_until;   // day confinement lifted per type (T+1 if never lifted)

  // ---- Change in notification-to-cull delay (13 Jan) ----
  int<lower=1> t_cap_change;

  // ---- Wild-bird pressure alpha_t ----
  // 1 before t_eps_ramp_start, exponential decay to t_eps_ramp_end, 0 after
  int<lower=1> t_eps_ramp_start;
  int<lower=1> t_eps_ramp_end;

  // ---- Initial states ----
  array[F] int<lower=1, upper=5> state_t0; // state at t=0: 1=S,2=E,3=I,4=N,5=R

  // ---- E duration: NegBin(sigma_r, sigma_p) ----
  real<lower=0> sigma_r;
  real<lower=0> sigma_p;

  // ---- Maximum E and I durations considered (days) ----
  int<lower=1> max_E_duration;

  int<lower=1> max_I_duration;
}

transformed data {
  // Likelihood cases:
  //   1: notified
  //   2: preventively culled in I
  //   3: preventively culled in E
  //   4: preventively culled in S or never infected

  array[F] int t_I_lo;   // lower bound on day farm entered I (0 if not applicable)
  array[F] int t_I_hi;   // upper bound on day farm entered I (0 if not applicable)
  array[F] int s_min;    // earliest candidate infection day
  array[F] int s_max;    // latest  candidate infection day

  // first_active[i]: first day farm i is active
  array[F] int first_active;
  for (i in 1:F) {
    first_active[i] = T + 1;
    for (t in 1:T) {
      if (active[i, t] > 0) { first_active[i] = t; break; }
    }
  }

  for (i in 1:F) {
    // ---- Case 1: notified farm ----
    if (t_notification[i] > 0) {
      t_I_hi[i] = t_notification[i] - 1;
      t_I_lo[i] = (t_suspicion[i] > 0) ? t_suspicion[i] : first_active[i];
      if (t_I_lo[i] >= t_notification[i])
        t_I_lo[i] = t_notification[i] - 1;
      // Cap: farm entered I at most max_I_duration days before notification.
      t_I_lo[i] = max(t_I_lo[i], t_notification[i] - max_I_duration);

    // ---- Case 2: prev culled in I ----
    // t_I entered before cull; suspicion (if present) gives lower bound
    } else if (t_prev_cull[i] > 0 && prev_cull_state[i] == 3) {
      t_I_hi[i] = t_prev_cull[i] - 1;
      t_I_lo[i] = (t_suspicion[i] > 0) ? t_suspicion[i] : first_active[i];
      if (t_I_lo[i] >= t_prev_cull[i])
        t_I_lo[i] = t_prev_cull[i] - 1;
      // Cap: farm entered I at most max_I_duration days before preventive cull.
      t_I_lo[i] = max(t_I_lo[i], t_prev_cull[i] - max_I_duration);

    // ---- Case 3: prev culled in E ----
    // Never entered I
    } else if (t_prev_cull[i] > 0 && prev_cull_state[i] == 2) {
      t_I_lo[i] = T + 1;   // sentinel: t_I never observed
      t_I_hi[i] = T + 1;

    // ---- Case 4: S at cull or never infected ----
    } else {
      t_I_lo[i] = 0; t_I_hi[i] = 0;
    }

    // ---- Compute s_min / s_max for cases 1, 2, 3 ----
    // Upper limit on infection day:
    //   cases 1 & 2: must be before t_I_lo (farm was already in I by then)
    //   case 3:      must be before t_prev_cull (farm was in E at cull)
    int s_ceil;
    if (t_notification[i] > 0 || prev_cull_state[i] == 3)
      s_ceil = t_I_lo[i] - 1;
    else if (prev_cull_state[i] == 2)
      s_ceil = t_prev_cull[i] - 1;
    else
      s_ceil = 0;

    if (s_ceil > 0) {
      // Lower limit: lookback from the latest plausible t_I entry
      int t_I_ref = (t_notification[i] > 0) ? t_I_hi[i]
                    : (prev_cull_state[i] == 3) ? t_I_hi[i]
                    : t_prev_cull[i] - 1;   // case 3: E duration ends at cull
      int raw_s_min = max(1, t_I_ref - max_E_duration);

      // Walk forward to find first active day in [raw_s_min, s_ceil]
      s_min[i] = 0;
      for (t in raw_s_min:s_ceil) {
        if (active[i, t] > 0) { s_min[i] = t; break; }
      }
      // Walk back from s_ceil to find last active day
      s_max[i] = 0;
      for (idx in 1:s_ceil) {
        int t = s_ceil - idx + 1;
        if (active[i, t] > 0) { s_max[i] = t; break; }
      }
      if (s_min[i] == 0) s_min[i] = s_max[i];
      if (s_max[i] < s_min[i]) s_max[i] = s_min[i];
    } else {
      s_min[i] = 0; s_max[i] = 0;
    }
  }

  // Confinement indicator [k, t]
  array[K, T] int confined;
  for (k in 1:K)
    for (t in 1:T)
      confined[k, t] = (t >= confined_from[k] && t < confined_until[k]) ? 1 : 0;

  // alpha_t: scaling of wild-bird pressure (data only)
  real eps_decay_rate = 5.0 / (t_eps_ramp_end - t_eps_ramp_start);
  vector[T] alpha_t;
  for (t in 1:T) {
    if (t < t_eps_ramp_start)
      alpha_t[t] = 1.0;
    else if (t < t_eps_ramp_end)
      alpha_t[t] = exp(-eps_decay_rate * (t - t_eps_ramp_start));
    else
      alpha_t[t] = 0.0;
  }

  // ---- Precomputations ----

  // 1 for premises with evidence of infection (cases 1-3)
  array[F] int has_infection;
  for (i in 1:F)
    has_infection[i] = (t_notification[i] > 0
                        || (t_prev_cull[i] > 0 && prev_cull_state[i] >= 2)) ? 1 : 0;

  // Potential sources: premises that are ever infectious or notified
  array[F] int inf_farms_buf;   // oversized buffer; first n_inf entries are valid
  int n_inf = 0;
  for (j in 1:F) {
    if (t_notification[j] > 0 || prev_cull_state[j] == 3
        || state_t0[j] == 3    || state_t0[j] == 4) {
      n_inf += 1;
      inf_farms_buf[n_inf] = j;
    }
  }

  // Receptors: premises in cases 1-3
  array[F] int recpt_buf;   // first n_recpt entries are valid
  int n_recpt = 0;
  for (i in 1:F)
    if (has_infection[i] == 1) {
      n_recpt += 1;
      recpt_buf[n_recpt] = i;
    }

  // Premises index -> receptor slot (0 for case 4)
  array[F] int recpt_idx;
  for (i in 1:F) recpt_idx[i] = 0;
  for (r in 1:n_recpt) recpt_idx[recpt_buf[r]] = r;

  // Width of the infection-window weight vector
  int W = max_E_duration;
}

parameters {
  real<lower=0>          gamma0;       // overall transmissibility
  real<lower=0, upper=1> gamma1;       // relative infectiousness: N vs I
  real<lower=0>          epsilon_hrz;
  real<lower=0>          epsilon_non_hrz;

  // Production-type weights (xi[1] = zeta[1] = 1); flock-size exponents fixed at 0.5
  vector<lower=0>[K-1]   xi_raw;
  vector<lower=0>[K-1]   zeta_raw;

  // Confinement multiplier on susceptibility and infectiousness
  real<lower=0, upper=1> tau;

  // Spatial kernel decay rate (km^-1)
  real<lower=0.01>       delta;

  // Mean infectious period before notification (days)
  real<lower=1, upper=7>  mu;
  real<lower=1, upper=60>  rho_n_1;     // mean notification-to-cull delay before 13 Jan
  real<lower=1, upper=20>  rho_n_2;     // mean notification-to-cull delay from 13 Jan
  real<lower=0, upper=10>  rho_cull;    // mean culling duration (days); lognormal
  real<lower=0.01, upper=10>  rho_cull_sd;  // sd of culling duration
}

transformed parameters {
  vector<lower=0>[K] xi;
  vector<lower=0>[K] zeta;
  xi[1]   = 1.0;
  zeta[1] = 1.0;
  for (k in 2:K) {
    xi[k]   = xi_raw[k-1];
    zeta[k] = zeta_raw[k-1];
  }

  // Wild-bird survival integral for case-4 premises
  vector[F] eps_active_sum;
  for (i in 1:F) {
    eps_active_sum[i] = 0.0;
    if (has_infection[i] == 0) {
      int tc = t_prev_cull[i];
      int ki = prod_type[i];
      for (t in 1:T) {
        if (active[i, t] > 0 && (tc == 0 || t < tc)) {
          real tau_eps_it = confined[ki, t] ? tau : 1.0;
          eps_active_sum[i] += alpha_t[t] * pow(active[i, t], 0.5) * tau_eps_it;
        }
      }
    }
  }

  // ---- Marginalised infection likelihood (survival + window weights) --------
  array[F] vector[W] log_w;
  for (i in 1:F) log_w[i] = rep_vector(negative_infinity(), W);

  // log_surv: case 4 initialised with the wild-bird term (farm-to-farm
  // pressure subtracted below); cases 1-3 accumulated over time
  vector[F] log_surv;
  for (i in 1:F) {
    if (has_infection[i] == 0) {
      int ki = prod_type[i];
      log_surv[i] = -(in_hrz[i] ? epsilon_hrz : epsilon_non_hrz) * xi[ki] * eps_active_sum[i];
    } else {
      log_surv[i] = 0.0;
    }
  }

  for (t in 1:T) {
    // Pressure on receptors (cases 1-3) from infectious and notified premises
    vector[n_recpt] pressure_I = rep_vector(0.0, n_recpt);
    vector[n_recpt] pressure_N = rep_vector(0.0, n_recpt);

    for (idx in 1:n_inf) {
      int j = inf_farms_buf[idx];
      if (active[j, t] <= 0) continue;

      int t_exit_I_j = (t_notification[j] > 0) ? t_notification[j]
                       : (prev_cull_state[j] == 3) ? t_prev_cull[j]
                       : T + 1;
      int in_I_j = (t_I_lo[j] > 0 && t_I_lo[j] <= T
                    && t >= t_I_lo[j] && t < t_exit_I_j) ? 1 : 0;
      // In N until culling starts
      int in_N_j = (t_notification[j] > 0
                    && t >= t_notification[j]
                    && (t_cull_start[j] == 0 || t < t_cull_start[j])) ? 1 : 0;
      // Premises notified before the start of the fitting period
      if (state_t0[j] == 4 && t_notification[j] == 0) {
        if (t_cull_start[j] == 0 || t < t_cull_start[j]) in_N_j = 1;
      }
      if (t == 1 && state_t0[j] == 3) in_I_j = 1;

      if (in_I_j == 1 || in_N_j == 1) {
        int kj = prod_type[j];
        real tau_j = confined[kj, t] ? tau : 1.0;
        real inf_j = gamma0 * zeta[kj] * pow(active[j, t], 0.5) * tau_j;

        for (p in j_ptr[j]:(j_ptr[j + 1] - 1)) {
          int i = nbr_i[p];
          if (active[i, t] <= 0) continue;
          int ki = prod_type[i];
          real tau_i = confined[ki, t] ? tau : 1.0;
          real kern_ij = exp(-delta * nbr_dist_km[p]);
          real beta_ij = inf_j * xi[ki] * pow(active[i, t], 0.5) * tau_i * kern_ij;

          if (has_infection[i] == 0) {
            // Case 4: update survival directly
            int tc_i = t_prev_cull[i];
            if (tc_i == 0 || t < tc_i)
              log_surv[i] -= (in_I_j == 1) ? beta_ij : gamma1 * beta_ij;
          } else {
            // Cases 1-3
            int ri = recpt_idx[i];
            if (in_I_j == 1) pressure_I[ri] += beta_ij;
            else              pressure_N[ri] += beta_ij;
          }
        }
      }
    }

    // Cases 1-3: survival and infection-day weights
    for (r in 1:n_recpt) {
      int i = recpt_buf[r];
      if (active[i, t] <= 0) continue;

      int t_cull_i = (t_notification[i] > 0) ? t_cull_start[i] : t_prev_cull[i];

      int still_S = (t < s_min[i])
                     && (t_cull_i == 0 || t < t_cull_i)
                     && state_t0[i] == 1
                     ? 1 : 0;
      int in_window = (s_min[i] > 0
                       && t >= s_min[i] && t <= s_max[i]
                       && (t_cull_i == 0 || t < t_cull_i)) ? 1 : 0;

      if (still_S == 1 || in_window == 1) {
        int ki = prod_type[i];
        real tau_eps_i  = confined[ki, t] ? tau : 1.0;
        real suscept_i  = xi[ki] * pow(active[i, t], 0.5) * tau_eps_i;
        real eps_i      = (in_hrz[i] ? epsilon_hrz : epsilon_non_hrz) * alpha_t[t] * suscept_i;
        real lam_i = eps_i + pressure_I[r] + gamma1 * pressure_N[r];
        log_surv[i] += -lam_i;

        if (in_window == 1) {
          int w = t - s_min[i] + 1;
          if (w >= 1 && w <= W) {
            // Survival to t - 1 times probability of infection on day t
            real log_surv_before = log_surv[i] + lam_i;
            real log_hit = (lam_i > 0) ? log1m_exp(-lam_i) : negative_infinity();
            if (t_I_lo[i] <= T) {
              int n_tI = t_I_hi[i] - t_I_lo[i] + 1;
              vector[n_tI] tI_lp;
              for (q in 1:n_tI) {
                int tI  = t_I_lo[i] + q - 1;
                int e_d = tI - t;
                real log_e = (e_d >= 0)
                  ? log_negbin_duration(e_d, sigma_r, sigma_p)
                  : negative_infinity();
                int i_lb = (t_notification[i] > 0)
                             ? t_notification[i] - tI
                             : t_prev_cull[i] - tI;
                real log_i_surv = (i_lb >= 1)
                  ? (i_lb - 1) * log1m(1.0 / mu)
                  : 0.0;
                tI_lp[q] = log_e + log_i_surv;
              }
              log_w[i][w] = log_surv_before + log_hit + log_sum_exp(tI_lp);
            } else {
              int e_cens = t_prev_cull[i] - t;
              real log_e_surv = negative_infinity();
              for (d in e_cens:max_E_duration * 3) {
                log_e_surv = log_sum_exp(log_e_surv,
                               log_negbin_duration(d, sigma_r, sigma_p));
              }
              log_w[i][w] = log_surv_before + log_hit + log_e_surv;
            }
          }
        }
      }
    }
  }
}

model {
  // ---- Priors ----
  gamma0          ~ gamma(3, 136);           // mean 0.022, sd 0.013
  gamma1          ~ beta(7, 4);              // mean 0.64, sd 0.14
  epsilon_hrz     ~ lognormal(log(0.0009), 0.6);
  epsilon_non_hrz ~ lognormal(log(0.00018), 0.6);
  // xi[2..5], zeta[2..5]: broiler stage two, layer, conventional duck, organic duck
  xi_raw[1]       ~ lognormal(log(1.5), 0.65);  // broiler stage two
  xi_raw[2]       ~ lognormal(log(0.78), 0.65); // layer
  xi_raw[3]       ~ lognormal(log(0.13), 0.6);  // conventional duck
  xi_raw[4]       ~ lognormal(log(0.4), 0.65);  // organic duck
  zeta_raw[1]     ~ lognormal(log(0.7), 0.65);  // broiler stage two
  zeta_raw[2]     ~ lognormal(log(0.38), 0.7);  // layer
  zeta_raw[3]     ~ lognormal(log(0.32), 0.7);  // conventional duck
  zeta_raw[4]     ~ lognormal(log(0.4), 0.7);   // organic duck
  tau             ~ beta(8, 2);               // mean 0.8, sd 0.12
  delta           ~ lognormal(log(0.17), 0.3);
  mu              ~ gamma(16, 4);            // mean 4, sd 1
  rho_n_1         ~ gamma(17, 5.9);          // mean 2.9, sd 0.7
  rho_n_2         ~ gamma(36, 15);           // mean 2.4, sd 0.4
  rho_cull        ~ gamma(2, 1);
  rho_cull_sd     ~ normal(0, 2);

  // ---- Likelihood ----

  // Part 1: notification-to-cull and culling durations
  for (i in 1:F) {
    // Notification-to-cull delay for reactively culled premises (delays of
    // 0 days are uninformative)
    if (t_notification[i] > 0 && t_cull_start[i] > 0 && t_prev_cull[i] == 0) {
      int dur_N = t_cull_start[i] - t_notification[i];
      if (dur_N >= 1) {
        real rho_n_i = rho_n_step(t_notification[i],
                                  t_cap_change, rho_n_1, rho_n_2);
        target += log_geom(dur_N, 1.0 / rho_n_i);
      }
    }
    // Culling duration
    if (t_cull_complete[i] > 0 && t_cull_start[i] > 0) {
      int dur_cull = t_cull_complete[i] - t_cull_start[i];
      if (dur_cull >= 1)
        target += lognormal_lpdf(dur_cull | log(rho_cull), rho_cull_sd);
    }
  }

  // Part 2: infection likelihood (fixed penalty when no infection day is
  // feasible, which avoids NaN gradients from log_sum_exp of all -Inf)
  for (i in 1:F) {
    if (t_notification[i] > 0
        || (t_prev_cull[i] > 0 && prev_cull_state[i] >= 2)) {
      real max_w = max(log_w[i]);
      if (!is_inf(max_w))
        target += log_sum_exp(log_w[i]);
      else
        target += -100.0;
    } else {
      target += log_surv[i];
    }
  }
}

generated quantities {
  // Per-premises log-likelihood (for LOO)
  vector[F] log_lik;
  for (i in 1:F) {
    if (t_notification[i] > 0
        || (t_prev_cull[i] > 0 && prev_cull_state[i] >= 2)) {
      real max_w = max(log_w[i]);
      log_lik[i] = is_inf(max_w) ? -100.0 : log_sum_exp(log_w[i]);
    } else {
      log_lik[i] = log_surv[i];
    }
  }
}
