%% creation of the pde

% create geometry
gd = [3,4,0,1,1,0,1,1,0,0]';
gm = [0,1,0,1]';
ns = [82,49]';
sf = 'R1';
g = decsg(gd,sf,ns);

% create the pde with domain
model = createpde();
geometryFromEdges(model,g);

% create the boundary conditions
applyBoundaryCondition(model,'neumann','Edge',[1,2,3],'g',0);
applyBoundaryCondition(model,'dirichlet','Edge',4,'u',0);

% source term
sfun = @(location,state) -location.y.*exp(-state.time);

% create coefficients
f = sfun;
a = 0;
c = 4;
m = 1;
specifyCoefficients(model,'m',m,'d',0,'c',c,'a',a,'f',f);

% initial condition function
ifun1 = @(location) 2*(location.x-1).^2.*cos(pi*location.y).*sin(0.5*pi*location.x);
ifun2 = @(location) -2*(asinh(0.5*location.x.*location.y) + erf(0.5*location.y))./(sqrt(location.x.^2+location.y.^2)+1);

% create initial condition
u0 = ifun1;
ut0 = ifun2;
setInitialConditions(model,u0,ut0);

%% solve the pde

% generate the mesh
generateMesh(model);

% indicate the time and spacing
n = 151;
tlist = linspace(0,5,n);

% solve the pde
model.SolverOptions.ReportStatistics ='on';
results = solvepde(model,tlist);

%% plot the solution

u = results.NodalSolution;

% plotting
f = figure;
umax = max(max(u));
umin = min(min(u));
for i = 1:n
    pdeplot(model,'XYData',u(:,i),'ZData',u(:,i),'ZStyle','continuous','Mesh','off');
    axis([0 1 0 1 umin umax]); 
    caxis([umin umax]);
    xlabel x
    ylabel y
    zlabel u
    M(i) = getframe;
end

%% interpolate the solution

% parameterize the circle
xint = 0.35*cos(linspace(0,2*pi,n)) + 0.5;
yint = 0.35*sin(linspace(0,2*pi,n)) + 0.5;

% interpolate
uintrp = interpolateSolution(results,xint,yint,1:n);

% extract the data points
uq = diag(uintrp);

% export the interpolation with data
rawdata = [xint',yint',uq,tlist'];
writematrix(rawdata,'data.csv')

%% plot the interpolation

% plot the interpolation
intp = figure;
plot3(xint,yint,uq','LineWidth',1);
exportgraphics(intp,'interpolated.pdf','ContentType','vector')

%% data normalization

unorm = 50 + (uq - min(uq))*(2950)/(max(uq)-min(uq));

datanorm = [xint',yint',unorm];

%% generate the audio

A = datanorm(:,1);       % amplitude
F = datanorm(:,3);       % frequency
s = 44100;              % sample rate
b = 24;                 % bits
t = linspace(0,1,s);  % time

% probability
p = datanorm(:,2);
P = [1-p, p];
val = [0 1];

% generate sound with probability and built-in fade out
for i = 1:151
    R = randsample(val,1,true,P(i,:));
    if R == 1
        R = F(i);
    end
    X(:,2*i-1) = (1 - t.^3)*A(i).*sin(2*pi*R*t);
    X(:,2*i) = 0;
end

% reshape matrix to a vector
x = reshape(X,1,[]);

% play/save the sound
% sound(x,s,b)
audiowrite('generated.wav',x,s,'BitsPerSample',24)