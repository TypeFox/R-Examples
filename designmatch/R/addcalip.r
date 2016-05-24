#! From Paul's Design of Obs. Studies (p. 251): "...adds a penalty to a distance matrix dmat for violations of the calip_size." "The default calip_size width is 0.2*sd(p). The magnitude of the penalty is penalty multiplied by the magnitude of the violation, where penalty is set to 1000 by default."
.addcalip = function (dist_mat, t_ind, calip_cov, calip_size, calip_penalty) {
    sdp = sd(calip_cov)
    penalty_mat = abs(outer(calip_cov[t_ind == 1], calip_cov[t_ind == 0], "-"))
    penalty_mat = (penalty_mat>(calip_size*sdp))*calip_penalty
    dist_mat = dist_mat+penalty_mat
    dist_mat
}