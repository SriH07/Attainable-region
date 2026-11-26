
function van_de_vusse_AR_final()

clc;
close all;
tStart = tic;

set(0,'DefaultFigureVisible','on');

%% ================= PARAMETERS =================

base = 39;
k_ref = base*[1 2 3 4];

T_ref = 350;
T = 350;

Ea = [40e3 50e3 60e3 45e3];
Rg = 8.314;

k = k_ref .* exp(-Ea/Rg .* (1/T - 1/T_ref));

k1 = k(1);
k2 = k(2);
k3 = k(3);

m4 = 2;
k4e = m4*k(4);

CA0 = 1;
CB0 = 0;

rhs = @(t,C) [-k1*C(1)+k2*C(2)-k4e*C(1)^2;
                k1*C(1)-(k2+k3)*C(2)];

cstr = @(tau,CAin,CBin) ...
    cstr_general(tau,k1,k2,k3,k4e,CAin,CBin);

%% ================= SETTINGS =================

tauLim = [1e-7 1e2];
tauEnd = 1;

Nb = 250;
Ng = 15;
nTau = 60;
nT = 40;

maxIter = 20;

tol.area = 1e-5;
tol.ext = 1e-4;
tol.clos = 1e-4;
tol.BC = 1e-3;

hmax = 1e-3;

pfrB = @(X0,tt) pfr_batch(k1,k2,k3,k4e,X0,tt,hmax);

%% ================= INITIAL SEED =================

fprintf('\nBuilding initial reactor loci...\n');

tgFine = [0 logspace(-6,log10(tauEnd),150)];

Cf = pfrB([CA0 CB0],tgFine);

nSeed = 300;

tauSeed = logspace(log10(tauLim(1)), ...
                   log10(tauLim(2)),nSeed);

[CAs,CBs] = cstr(tauSeed,CA0,CB0);

pts = [Cf; CAs(:) CBs(:); CA0 CB0; 0 0];

P = hull_pts(pts);

% Save initial mixing region
P0 = P;

% Placeholder for region after first iteration
Pafter1 = P0;

%% ================= AR CONSTRUCTION =================

fprintf('\nStarting attainable-region construction...\n');

tau_c = logspace(log10(tauLim(1)), ...
                 log10(tauLim(2)),nTau);

tg = [0 logspace(-5,log10(tauEnd),nT)];

areaOld = polyarea(P(:,1),P(:,2));

hist = zeros(maxIter,3);
converged = false;

for iter = 1:maxIter

    B = resample_boundary(P,Nb);
    G = interior_grid(P,Ng);

    S = [B;G];

    % PFR outputs
    newP = pfrB(S,tg);

    % CSTR outputs
    [CAc,CBc] = cstr(tau_c,S(:,1),S(:,2));

    newC = [CAc(:) CBc(:)];

    % New reactor-output hull
    new = hull_pts([newP;newC]);

    % Extension beyond old region
    ext = max_violation(new,P);

    % Mixing closure
    P = hull_pts([P(1:end-1,:);new]);

    % Save region after first iteration
    if iter == 1
        Pafter1 = P;
    end

    areaNew = polyarea(P(:,1),P(:,2));

    relArea = abs(areaNew-areaOld)/max(areaNew,eps);

    hist(iter,:) = [areaNew relArea ext];

    fprintf(['Iteration %2d: Area = %.8f, ' ...
             'Relative change = %.2e, Extension = %.2e\n'], ...
             iter,areaNew,relArea,ext);

    areaOld = areaNew;

    if relArea < tol.area && ext < tol.ext

        converged = true;
        hist = hist(1:iter,:);
        break;

    end

end

if ~converged

    hist = hist(1:iter,:);

    warning('Construction did not meet both convergence tolerances.');

end

%% ================= INITIAL VS UPDATED AREA =================

area0 = polyarea(P0(:,1),P0(:,2));
area1 = polyarea(Pafter1(:,1),Pafter1(:,2));

fprintf('\n============================================\n');
fprintf('       INITIAL VS UPDATED MIXING REGION\n');
fprintf('============================================\n');

fprintf('Initial area        = %.10f\n',area0);
fprintf('After first update  = %.10f\n',area1);
fprintf('Area increase       = %.10e\n',area1-area0);

if area0 > 0
    fprintf('Percentage increase = %.6f %%\n', ...
        100*(area1-area0)/area0);
end

fprintf('Initial vertices    = %d\n',size(P0,1)-1);
fprintf('Updated vertices    = %d\n',size(Pafter1,1)-1);

%% ================= CLOSURE VERIFICATION =================

fprintf('\nPerforming closure verification...\n');

rng(2);

nV = size(P,1)-1;

Vx = P(1:end-1,1);
Vy = P(1:end-1,2);

nRandom = 500;

idx = randi(nV,nRandom,3);

W = -log(rand(nRandom,3));
W = W./sum(W,2);

Xr = [sum(W.*Vx(idx),2), ...
      sum(W.*Vy(idx),2)];

Xc = [interior_grid(P,25);
      resample_boundary(P,300);
      Xr;
      CA0 CB0;
      0 0];

Xp = [interior_grid(P,12);
      resample_boundary(P,100);
      Xr(1:200,:)];

% CSTR verification
tauV = logspace(log10(tauLim(1)), ...
                log10(tauLim(2)),301);

[CAv,CBv] = cstr(tauV,Xc(:,1),Xc(:,2));

Qc = [CAv(:) CBv(:)];

extC = max_violation(Qc,P);

% PFR verification
tVerify = [0 logspace(-5,log10(tauEnd),100)];

Qp = pfrB(Xp,tVerify);

extP = max_violation(Qp,P);

fprintf('CSTR closure violation = %.3e\n',extC);
fprintf('PFR closure violation  = %.3e\n',extP);

%% ================= BOUNDARY DIAGNOSTIC =================

fprintf('\nCalculating boundary-vector diagnostic...\n');

worst = -inf;

for e = 1:size(P,1)-1

    d = P(e+1,:)-P(e,:);

    if norm(d)<1e-12
        continue;
    end

    nrm = [d(2) -d(1)]/norm(d);

    if abs(nrm(1))>0.999 && P(e,1)<1e-9
        continue;
    end

    if abs(nrm(2))>0.999 && P(e,2)<1e-9
        continue;
    end

    for t = linspace(0.05,0.95,8)

        x = P(e,:)+t*d;

        f = rhs(0,x(:));

        outward = dot(f,nrm)/max(norm(f),1e-12);

        worst = max(worst,outward);

    end

end

if ~isfinite(worst)
    worst = inf;
end

%% ================= FINAL RESULTS =================

area = polyarea(P(:,1),P(:,2));

[CBmax,im] = max(P(:,2));
CAatCBmax = P(im,1);

fprintf('\n');
fprintf('============================================\n');
fprintf('       FINAL ATTAINABLE REGION RESULTS\n');
fprintf('============================================\n');

fprintf('Temperature           : %.2f K\n',T);

fprintf('k1                    : %.4f\n',k1);
fprintf('k2                    : %.4f\n',k2);
fprintf('k3                    : %.4f\n',k3);
fprintf('k4                    : %.4f\n',k(4));

fprintf('Region area           : %.8f (mol/L)^2\n',area);

fprintf('Maximum CB            : %.6f mol/L\n',CBmax);
fprintf('CA at maximum CB      : %.6f mol/L\n',CAatCBmax);

fprintf('Boundary vertices     : %d\n',size(P,1)-1);

fprintf('Construction converged: %d\n',converged);

fprintf('CSTR closure violation: %.3e\n',extC);
fprintf('PFR closure violation : %.3e\n',extP);

fprintf('Boundary diagnostic   : %.3e\n',worst);

ok = converged && ...
     max(extC,extP)<tol.clos && ...
     worst<tol.BC;

if ok
    fprintf('\nVERDICT: All specified numerical checks passed.\n');
else
    fprintf('\nVERDICT: One or more numerical checks failed.\n');
end

fprintf('Total execution time: %.2f seconds\n',toc(tStart));

%% =========================================================
% CALCULATE FOUR-SPECIES CONCENTRATION PROFILES
% These calculations do not modify the AR construction.
%% =========================================================

fprintf('\nCalculating four-species concentration profiles...\n');

% PFR profiles
tau_pfr = linspace(0,tauEnd,600);

[~,Cp] = ode45( ...
    @(t,C) full_rhs(t,C,k1,k2,k3,k(4)), ...
    tau_pfr,[CA0;CB0;0;0]);

% CSTR profiles
tau_cstr_profile = logspace( ...
    log10(tauLim(1)),log10(tauLim(2)),500);

[CAc_profile,CBc_profile] = cstr( ...
    tau_cstr_profile,CA0,CB0);

% CSTR product balances
CCc_profile = tau_cstr_profile(:).*k3.*CBc_profile(:);

CDc_profile = tau_cstr_profile(:).*k(4).*CAc_profile(:).^2;

%% ================= COMMON PLOT LIMITS =================

xmax = 1.05*CA0;
ymax = max(1.2*CBmax,0.05);

%% ================= PLOT 1: PFR PROFILES =================

f1 = figure(1);
clf(f1);
set(f1,'Color','w','Visible','on');

ax1 = axes('Parent',f1);

plot(ax1,tau_pfr,Cp(:,1),'b-','LineWidth',2);
hold(ax1,'on');

plot(ax1,tau_pfr,Cp(:,2),'r-','LineWidth',2);
plot(ax1,tau_pfr,Cp(:,3),'g-','LineWidth',2);
plot(ax1,tau_pfr,Cp(:,4),'m-','LineWidth',2);

xlabel(ax1,'PFR residence time, \tau (s)');
ylabel(ax1,'Concentration (mol/L)');

title(ax1,'PFR Concentration Profiles');

legend(ax1,'C_A','C_B','C_C','C_D', ...
    'Location','best');

grid(ax1,'on');
box(ax1,'on');

%% ================= PLOT 2: CSTR PROFILES =================

f2 = figure(2);
clf(f2);
set(f2,'Color','w','Visible','on');

ax2 = axes('Parent',f2);

semilogx(ax2,tau_cstr_profile,CAc_profile,'b-','LineWidth',2);
hold(ax2,'on');

semilogx(ax2,tau_cstr_profile,CBc_profile,'r-','LineWidth',2);
semilogx(ax2,tau_cstr_profile,CCc_profile,'g-','LineWidth',2);
semilogx(ax2,tau_cstr_profile,CDc_profile,'m-','LineWidth',2);

xlabel(ax2,'CSTR residence time, \tau (s)');
ylabel(ax2,'Concentration (mol/L)');

title(ax2,'CSTR Concentration Profiles');

legend(ax2,'C_A','C_B','C_C','C_D', ...
    'Location','best');

grid(ax2,'on');
box(ax2,'on');

%% ================= PLOT 3: OPERATING TRAJECTORIES =================

f3 = figure(3);
clf(f3);
set(f3,'Color','w','Visible','on');

ax3 = axes('Parent',f3);

plot(ax3,Cf(:,1),Cf(:,2),'b-','LineWidth',1.8);
hold(ax3,'on');

plot(ax3,CAs(:),CBs(:),'r--','LineWidth',1.8);

plot(ax3,CA0,CB0,'ko','MarkerFaceColor','k');

% Direction arrows along increasing residence time
add_curve_arrows(ax3,Cf,'b',8);

CSTRpoints = [CAs(:) CBs(:)];

add_curve_arrows(ax3,CSTRpoints,'r',8);

xlabel(ax3,'C_A (mol/L)');
ylabel(ax3,'C_B (mol/L)');

title(ax3,'PFR and CSTR Operating Trajectories');

legend(ax3,'PFR trajectory','CSTR trajectory', ...
    'Feed','Location','best');

xlim(ax3,[0 xmax]);
ylim(ax3,[0 ymax]);

grid(ax3,'on');
box(ax3,'on');

%% ================= PLOT 4: INITIAL MIXING REGION =================

f4 = figure(4);
clf(f4);
set(f4,'Color','w','Visible','on');

ax4 = axes('Parent',f4);

fill(ax4,P0(:,1),P0(:,2),[0.85 0.92 1], ...
    'EdgeColor','b','LineWidth',1.8);

hold(ax4,'on');

plot(ax4,Cf(:,1),Cf(:,2),'b-','LineWidth',1.2);
plot(ax4,CAs(:),CBs(:),'r--','LineWidth',1.2);

plot(ax4,CA0,CB0,'ko','MarkerFaceColor','k');

xlabel(ax4,'C_A (mol/L)');
ylabel(ax4,'C_B (mol/L)');

title(ax4,'Initial Mixing Region');

legend(ax4,'Initial mixing region','PFR locus', ...
    'CSTR locus','Feed','Location','best');

xlim(ax4,[0 xmax]);
ylim(ax4,[0 ymax]);

grid(ax4,'on');
box(ax4,'on');

%% ================= PLOT 5: UPDATED MIXING REGION =================

f5 = figure(5);
clf(f5);
set(f5,'Color','w','Visible','on');

ax5 = axes('Parent',f5);

fill(ax5,Pafter1(:,1),Pafter1(:,2),[0.85 0.92 1], ...
    'EdgeColor','b','LineWidth',1.8);

hold(ax5,'on');

plot(ax5,P0(:,1),P0(:,2),'r--','LineWidth',1.5);

xlabel(ax5,'C_A (mol/L)');
ylabel(ax5,'C_B (mol/L)');

title(ax5,'Updated Mixing Region After First Iteration');

legend(ax5,'Updated mixing region','Initial mixing region', ...
    'Location','best');

xlim(ax5,[0 xmax]);
ylim(ax5,[0 ymax]);

grid(ax5,'on');
box(ax5,'on');

%% ================= PLOT 6: FINAL ATTAINABLE REGION =================

f6 = figure(6);
clf(f6);
set(f6,'Color','w','Visible','on');

ax6 = axes('Parent',f6);

plot(ax6,P(:,1),P(:,2),'k-','LineWidth',2);
hold(ax6,'on');

plot(ax6,Cf(:,1),Cf(:,2),'b-','LineWidth',1.8);

plot(ax6,CAs(:),CBs(:),'r--','LineWidth',1.8);

plot(ax6,CA0,CB0,'ko','MarkerFaceColor','k');

plot(ax6,CAatCBmax,CBmax,'gp', ...
    'MarkerFaceColor','g','MarkerSize',10);

% Direction arrows on reactor trajectories
add_curve_arrows(ax6,Cf,'b',8);

add_curve_arrows(ax6,CSTRpoints,'r',8);

xlabel(ax6,'C_A (mol/L)');
ylabel(ax6,'C_B (mol/L)');

title(ax6,'Final Van de Vusse Attainable Region');

legend(ax6,'Attainable region','PFR locus','CSTR locus', ...
    'Feed','Maximum C_B','Location','best');

xlim(ax6,[0 xmax]);
ylim(ax6,[0 ymax]);

grid(ax6,'on');
box(ax6,'on');

%% ================= PLOT 7: CONVERGENCE =================

f7 = figure(7);
clf(f7);
set(f7,'Color','w','Visible','on');

ax7 = axes('Parent',f7);

iterations = 1:size(hist,1);

plot(ax7,iterations,hist(:,2),'bo-', ...
    'LineWidth',1.5,'MarkerFaceColor','b');

hold(ax7,'on');

plot(ax7,iterations,hist(:,3),'rs-', ...
    'LineWidth',1.5,'MarkerFaceColor','r');

xlabel(ax7,'Iteration');
ylabel(ax7,'Value');

title(ax7,'Attainable Region Convergence');

legend(ax7,'Relative area change','Boundary extension');

grid(ax7,'on');
box(ax7,'on');

drawnow;

fprintf('\nAll calculations and figures displayed.\n');

end


%% ============================================================
% FOUR-SPECIES PFR DIFFERENTIAL EQUATIONS
%% ============================================================

function dC = full_rhs(~,C,k1,k2,k3,k4)

CA = C(1);
CB = C(2);

dC = zeros(4,1);

dC(1) = -k1*CA+k2*CB-2*k4*CA^2;
dC(2) = k1*CA-(k2+k3)*CB;
dC(3) = k3*CB;
dC(4) = k4*CA^2;

end


%% ============================================================
%                     ORIGINAL LOCAL FUNCTIONS
%% ============================================================

function Q = pfr_batch(k1,k2,k3,k4e,X0,tt,hmax)

if isvector(X0)
    X0 = reshape(X0,1,2);
end

N = size(X0,1);
nt = numel(tt);

blk = 1000;

Q = zeros(N*nt,2);
pos = 0;

for ib = 1:blk:N

    X = X0(ib:min(ib+blk-1,N),:);

    n = size(X,1);

    Qb = zeros(n*nt,2);

    Qb(1:n,:) = X;

    for j = 2:nt

        dt = tt(j)-tt(j-1);

        ns = max(1,ceil(dt/hmax));

        h = dt/ns;

        for s = 1:ns

            K1 = rate(X,k1,k2,k3,k4e);

            K2 = rate(X+0.5*h*K1,k1,k2,k3,k4e);

            K3 = rate(X+0.5*h*K2,k1,k2,k3,k4e);

            K4 = rate(X+h*K3,k1,k2,k3,k4e);

            X = X+(h/6)*(K1+2*K2+2*K3+K4);

            if any(X(:)<-1e-8)
                error('RK4 generated significantly negative concentrations.');
            end

            X(X<0) = 0;

        end

        rows = (j-1)*n+(1:n);

        Qb(rows,:) = X;

    end

    m = n*nt;

    Q(pos+1:pos+m,:) = Qb;

    pos = pos+m;

end

Q = Q(1:pos,:);

end


function F = rate(X,k1,k2,k3,k4e)

CA = X(:,1);
CB = X(:,2);

dCA = -k1*CA+k2*CB-k4e*CA.^2;

dCB = k1*CA-(k2+k3)*CB;

F = [dCA dCB];

end


function [CA,CB] = cstr_general(tau,k1,k2,k3,k4e,CAin,CBin)

den = 1+(k2+k3)*tau;

alpha = k1*tau./den;

beta = CBin./den;

a = k4e*tau;

b = 1+tau.*(k1-k2*alpha);

c = CAin+tau*k2.*beta;

CA = 2*c./(b+sqrt(b.^2+4*a.*c));

CB = beta+alpha.*CA;

end


function H = hull_pts(Q)

if isempty(Q)
    H = zeros(0,2);
    return;
end

blk = 1e5;

H = zeros(0,2);

for i = 1:blk:size(Q,1)

    q = Q(i:min(i+blk-1,size(Q,1)),:);

    q = try_hull(q);

    H = [H;q];

end

H = try_hull(H);

end


function q = try_hull(q)

if size(q,1)>=3

    try
        K = convhull(q(:,1),q(:,2));
        q = q(K,:);
    catch
    end

end

end


function G = interior_grid(P,N)

[gx,gy] = meshgrid( ...
    linspace(min(P(:,1)),max(P(:,1)),N), ...
    linspace(min(P(:,2)),max(P(:,2)),N));

in = inpolygon(gx(:),gy(:),P(:,1),P(:,2));

G = [gx(in) gy(in)];

end


function B = resample_boundary(P,Nb)

seg = sqrt(sum(diff(P).^2,2));

s = [0;cumsum(seg)];

[s,iu] = unique(s);

Pu = P(iu,:);

if numel(s)<2 || s(end)<=0
    B = P(1:end-1,:);
    return;
end

sq = linspace(0,s(end),Nb+1);

sq(end) = [];

B = [interp1(s,Pu(:,1),sq(:)), ...
     interp1(s,Pu(:,2),sq(:))];

B = unique([B;P(1:end-1,:)],'rows');

end


function v = max_violation(Q,P)

Qh = hull_pts(Q);

if isempty(Qh)
    v = 0;
    return;
end

out = ~inpolygon(Qh(:,1),Qh(:,2),P(:,1),P(:,2));

v = 0;

if any(out)
    v = max(dist_to_polygon(Qh(out,:),P));
end

end


function d = dist_to_polygon(Q,P)

d = inf(size(Q,1),1);

for e = 1:size(P,1)-1

    A = P(e,:);
    B = P(e+1,:);

    AB = B-A;

    L2 = max(AB*AB',1e-24);

    t = ((Q(:,1)-A(1))*AB(1)+ ...
         (Q(:,2)-A(2))*AB(2))/L2;

    t = min(max(t,0),1);

    d = min(d,hypot( ...
        Q(:,1)-(A(1)+t*AB(1)), ...
        Q(:,2)-(A(2)+t*AB(2))));

end

end


%% ============================================================
% DIRECTION ARROWS ALONG OPERATING TRAJECTORIES
%% ============================================================

function add_curve_arrows(ax,X,color,nArrows)

X = X(all(isfinite(X),2),:);

if size(X,1)<5
    return;
end

idx = round(linspace(3,size(X,1)-2,nArrows));

dX = X(idx+1,:)-X(idx-1,:);

quiver(ax, ...
    X(idx,1),X(idx,2), ...
    dX(:,1),dX(:,2), ...
    0.5, ...
    'Color',color, ...
    'LineWidth',1.4, ...
    'MaxHeadSize',2.5);

end
