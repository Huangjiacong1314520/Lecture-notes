function out = skew_mu_lower_bound_Holland(M, k_f, k_v, opts)
if nargin < 4, opts = struct; end
if ~isfield(opts,'maxIter'), opts.maxIter = 200; end
if ~isfield(opts,'tol'), opts.tol = 1e-6; end
if ~isfield(opts,'verbose'), opts.verbose = 1; end

nf = sum(k_f); nv = sum(k_v); n = nf+nv;
[nr,nc] = size(M);
if nr~=n || nc~=n, error('M size %dx%d inconsistent with nf+nv=%d.',nr,nc,n); end

rng(1);
b = randn(n,1)+1j*randn(n,1); b = b/norm(b);
w = randn(n,1)+1j*randn(n,1); w = w/norm(w);
v_total = 1; flag = 0; M_scaled = M;

hist_v=[]; hist_g1=[]; hist_c1=[];

for iter=1:opts.maxIter
    % ---- RIGHT: b -> a ----
    g = M_scaled*b; g1=g(1:nf); g2=g(nf+1:end);
    n_g1 = norm(g1); n_g2 = norm(g2);
    if n_g1 >= 1-1e-9, flag=1; if opts.verbose, fprintf('[%3d] ||g1||>=1\n',iter); end; break; end
    v_star = (n_g2^2)/sqrt(max(1e-16,1-n_g1^2));
    S_star = blkdiag(eye(nf), v_star*eye(nv));
    a = S_star \ g; a = a/norm(a);

    % ---- LEFT: z -> w ----
    z = local_update_z(a,w,k_f,k_v,nf,nv);
    c = M_scaled' * z; c1=c(1:nf); c2=c(nf+1:end);
    n_c1 = norm(c1); n_c2 = norm(c2);
    if n_c1 >= 1-1e-9, flag=1; if opts.verbose, fprintf('[%3d] ||c1||>=1\n',iter); end; break; end
    v_sharp = (n_c2^2)/sqrt(max(1e-16,1-n_c1^2));
    S_sharp = blkdiag(eye(nf), v_sharp*eye(nv));
    w_new = S_sharp \ c; w_new = w_new/norm(w_new);

    % ---- combine & scale ----
    v_local = max(v_star, v_sharp);
    hist_v(end+1)=v_local; hist_g1(end+1)=n_g1; hist_c1(end+1)=n_c1;

    if v_local < 1+opts.tol
        if opts.verbose, fprintf('[%3d] v<1+tol -> STOP, lb=%.6g\n',iter,v_total); end
        break;
    end
    v_total = v_total * v_local;

    S_L = blkdiag(eye(nf), 1/sqrt(v_local)*eye(nv));
    S_R = blkdiag(eye(nf), sqrt(v_local)*eye(nv));
    M_scaled = S_L * M_scaled * S_R;

    % update right vector using block rules
    b = local_update_b(a, w_new, k_f, k_v, nf, nv); b = b/norm(b);
    w = w_new;

    if v_total>1e10, flag=1; break; end
end

out.mu_lb   = v_total;
out.v_total = v_total;
out.iter    = iter;
out.flag    = flag;
out.history = struct('v',hist_v,'g1',hist_g1,'c1',hist_c1);
end

% ---- helpers ----
function z = local_update_z(a,w,kf,kv,nf,nv)
z=zeros(nf+nv,1); idx=1;
if kf(1)>0
    a1=a(idx:idx+kf(1)-1); w1=w(idx:idx+kf(1)-1);
    z(idx:idx+kf(1)-1)=sign(real(a1.*conj(w1))).*w1; idx=idx+kf(1);
end
if kf(2)>0
    a2=a(idx:idx+kf(2)-1); w2=w(idx:idx+kf(2)-1);
    inner=a2'*w2; z(idx:idx+kf(2)-1)=(inner/abs(inner+eps))*w2; idx=idx+kf(2);
end
if kf(3)>0
    a3=a(idx:idx+kf(3)-1); w3=w(idx:idx+kf(3)-1);
    z(idx:idx+kf(3)-1)=(abs(w3)./max(abs(a3),1e-12)).*a3; idx=idx+kf(3);
end
for blk=1:3
    start=nf+sum(kv(1:blk-1))+1; stop=start+kv(blk)-1;
    if kv(blk)==0, continue; end
    av=a(start:stop); wv=w(start:stop);
    switch blk
        case 1
            z(start:stop)=sign(real(av.*conj(wv))).*wv;
        case 2
            inner=av'*wv; z(start:stop)=(inner/abs(inner+eps))*wv;
        case 3
            z(start:stop)=(abs(wv)./max(abs(av),1e-12)).*av;
    end
end
end

function b = local_update_b(a,w,kf,kv,nf,nv)
b=zeros(nf+nv,1); idx=1;
if kf(1)>0
    a1=a(idx:idx+kf(1)-1); w1=w(idx:idx+kf(1)-1);
    b(idx:idx+kf(1)-1)=sign(real(a1.*conj(w1))).*a1; idx=idx+kf(1);
end
if kf(2)>0
    a2=a(idx:idx+kf(2)-1); w2=w(idx:idx+kf(2)-1);
    inner=a2'*w2; b(idx:idx+kf(2)-1)=(inner/abs(inner+eps))*a2; idx=idx+kf(2);
end
if kf(3)>0
    a3=a(idx:idx+kf(3)-1); w3=w(idx:idx+kf(3)-1);
    b(idx:idx+kf(3)-1)=(abs(a3)./max(abs(w3),1e-12)).*w3; idx=idx+kf(3);
end
for blk=1:3
    start=nf+sum(kv(1:blk-1))+1; stop=start+kv(blk)-1;
    if kv(blk)==0, continue; end
    av=a(start:stop); wv=w(start:stop);
    switch blk
        case 1
            b(start:stop)=sign(real(av.*conj(wv))).*av;
        case 2
            inner=av'*wv; b(start:stop)=(inner/abs(inner+eps))*av;
        case 3
            b(start:stop)=(abs(av)./max(abs(wv),1e-12)).*wv;
    end
end
end
