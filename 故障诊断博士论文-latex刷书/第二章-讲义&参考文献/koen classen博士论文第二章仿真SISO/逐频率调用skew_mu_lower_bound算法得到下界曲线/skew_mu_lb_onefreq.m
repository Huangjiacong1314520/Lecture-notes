function out = skew_mu_lb_onefreq(Mjw, k_f, k_v, opts)
% 单频点：对 M(jw) 调用 Holland-style 幂迭代
% 需要你把“附：幂迭代实现”中的 skew_mu_lower_bound_Holland.m 放在路径上
out = skew_mu_lower_bound_Holland(Mjw, k_f, k_v, opts);
end
