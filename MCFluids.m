%% Initialize 
global n_grid t_max w_t vx vy dt nd dx nu filled calls WoS_eps;
n_grid = 32;
t_max = 100;
filled = 0;
calls = 0;
WoS_eps = 0.1;

%initial conditions for vorticity
w_t = NaN(n_grid, n_grid, t_max + 1);
vx =  NaN(n_grid, n_grid, t_max + 1);
vy =  NaN(n_grid, n_grid, t_max + 1);

% w_t(:,:, 1) = zeros(n_grid, n_grid);
% w_t([8:14], [8:14], 1) = 1;
% w_t([18:24], [18:24], 1) = -1;

% w_t(:,:, 1) = zeros(n_grid, n_grid);
% w_t([35:55], [35:55], 1) = 1;
% w_t([75:105], [75:105], 1) = -1;


for i=1:n_grid
    for j=1:n_grid
        if mod(i+j, 2) == 0 
            w_t(i,j,1) = 1;
        else
            w_t(i,j,1) = -1;
        end
    end
end


dt = 0.01;
nd = 10;
dx = 1;
nu = 100;
%% Compute vorticity
for t_step = 1:t_max
    for i = 0:(n_grid-1)
        for j = 0:(n_grid-1)
            vort(t_step*dt, dx*[i, j]);
        end
    end
end
%disp(w_t(:,:, 3));

%% Compute velocities from the saved vorticity array
for t = 0:t_max
    for i = 0:(n_grid-1)
        for j = 0:(n_grid-1)
            %We can use either Biot-Savart or WoS to compute the velocities
            v = compvelWoS(t*dt, [i, j]);
            %v = compvelBS(t*dt, [i, j]);
            vx(i+1, j+1, t+1) = v(1);
            vy(i+1, j+1, t+1) = v(2);
        end
    end
    disp(t);
end


%% Display vorticity animation
for i = 1:(t_max + 1)
    imagesc(w_t(:,:,i));
    title(["t = ", dt*i]);
    drawnow;
end

%%  Save vorticity animation to gif
% nImages = t_max+1;
% 
% fig = figure;
% for idx = 1:nImages
%     imagesc(w_t(:,:,idx), [-1, 1])
%     colorbar
%     title(["t = ", dt*idx])
%     drawnow
%     frame = getframe(fig);
%     im{idx} = frame2im(frame);
% end
% 
% filename = "t500_twoblob_128_nd10.gif"; % Specify the output file name
% for idx = 1:nImages
%     [A,map] = rgb2ind(im{idx},256);
%     if idx == 1
%         imwrite(A,map,filename,"gif","LoopCount",Inf,"DelayTime",0.03);
%     else
%         imwrite(A,map,filename,"gif","WriteMode","append","DelayTime",0.03);
%     end
%     disp(idx);
% end
% 
%% Display tracer particle animation
xs_all = zeros(floor((n_grid^2 - 1)/63), 2);

for j=0:15:(n_grid^2 - 1)
    xs_all(j+1, 1) = floor(j/n_grid);
    xs_all(j+1, 2) = mod(j, n_grid);
end


disp("made it");
for t_ind = 1:t_max
    for i = 1:15:(n_grid^2)
    %for i = 1:7
        xcur = [xs_all(i, 1); xs_all(i, 2)];
        x_nex = xcur + [vx(ind(xcur(1)), ind(xcur(2)), t_ind); vy(ind(xcur(1)), ind(xcur(2)), t_ind)];
        
        
        if (x_nex(1) >= 0) && (x_nex(1) <= (n_grid-1)*dx) && (x_nex(2) >= 0) && (x_nex(2) <= (n_grid-1)*dx )
        xs_all(i, :) = x_nex;
        end
    end
    imagesc(w_t(:,:,t_ind), [-1, 1])
    colorbar
    hold on
    scatter( xs_all(:,1), xs_all(:,2), "red", ".")
    axis([0, n_grid, 0, n_grid])
    title(["t = ", dt*t_ind])
    drawnow
    hold off
end

%%  Save tracer animation to gif

% nImages = t_max+1;
% 
% for j=0:15:(n_grid^2 - 1)
%     xs_all(j+1, 1) = floor(j/n_grid);
%     xs_all(j+1, 2) = mod(j, n_grid);
% end
% 
% fig = figure;
% for t_ind = 1:nImages
%     for i = 1:15:(n_grid^2)
%         xcur = [xs_all(i, 1); xs_all(i, 2)];
%         x_nex = xcur + [vx(ind(xcur(1)), ind(xcur(2)), t_ind); vy(ind(xcur(1)), ind(xcur(2)), t_ind)];
%         if (x_nex(1) >= 0) && (x_nex(1) <= (n_grid-1)*dx) && (x_nex(2) >= 0) && (x_nex(2) <= (n_grid-1)*dx )
%             xs_all(i, :) = x_nex;
%         end
%     end
%     imagesc(w_t(:,:,t_ind), [-1, 1])
%     colorbar
%     hold on
%     scatter( xs_all(:,1), xs_all(:,2), "red", ".")
%     axis([0, n_grid, 0, n_grid])
%     title(["t = ", dt*t_ind])
%     drawnow
%     frame = getframe(fig);
%     im{t_ind} = frame2im(frame);
%     hold off
% end
% 
% filename = "backtracksmall.gif"; % Specify the output file name
% for idx = 1:nImages
%     [A,map] = rgb2ind(im{idx},256);
%     if idx == 1
%         imwrite(A,map,filename,"gif","LoopCount",Inf,"DelayTime",0.03);
%     else
%         imwrite(A,map,filename,"gif","WriteMode","append","DelayTime",0.03);
%     end
%     disp(idx);
%  end


%%

%Biot-Savart velocity computer
function v_out = compvelBS(t, x)
global nd n_grid nu dt;
    sum = [0, 0];
    sig = n_grid/2; 
    for i = 1:nd
        y = normrnd(x, sig);
        g = (x - y)/(2*pi*(norm(x-y)^2));
        w = vort(t, [ ind(y(1)), ind(y(2)) ] );
        %w = w_t(ind(y(1)), ind(y(2)), indt(t));
        p = normpdf(y(1), x(1), sig)*normpdf(y(2), x(2), sig);

        sum = sum + [g(2)*w, - g(1)*w]/p;
    end
    v_out = sum/nd;
end

%WoS Velocity computer
function v = compvelWoS(t, x)
    grad = WoSgrad(t,x);
    v = [grad(2), -grad(1)];
end

%Backtracing, an alternate vorticity method presented in the paper, but not explored
%in my project
function w = backtrace(t, x)
global w_t dt calls;
    if isnan(w_t(ind(x(1)), ind(x(2)), indt(t)))
        calls = calls + 1;
        if mod(calls,1000) == 0
           disp(calls);
        end
        v = compvelBS(t - dt, x);
        w = backtrace(t - dt, x - v*dt);
        w_t(ind(x(1)), ind(x(2)), indt(t)) = w;
    else
        w = w_t(ind(x(1)), ind(x(2)), indt(t));
    end
end

% Get the vorticity at a point. Calls recursively back in time. Follows
% algorithm in my paper
function w = vort(t, x)
global w_t dt nd calls n_grid dx nu;
    if isnan(w_t(ind(x(1)), ind(x(2)), indt(t)))
        calls = calls + 1;
        if mod(calls,1000) == 0
            disp(calls);
        end
        
        sig = sqrt(2*nu*dt);
        
        w = 0;
        xa = x - dt*compvelWoS(t-dt, x);
        %xa = x - dt*compvelBS(t-dt, x);
        x_samp = zeros(2);
        for i = 1:nd
            %disp(xa);
            x_samp(1) = max([min([normrnd(xa(1), sig), (n_grid-1)*dx]), 0]);
            x_samp(2) = max([min([normrnd(xa(2), sig), (n_grid-1)*dx]), 0]);
            %Dw = stretch(t - dt, x_samp(i, :));
            w = w + (vort(t - dt, x_samp )); % + dt*Dw);
            
        end
        w = w/nd;
        w_t(ind(x(1)), ind(x(2)), indt(t)) = w;
    else
        if t == 0
            w = w_t(ind(x(1)),ind(x(2)), 1);
            return;
        end
        w = w_t(ind(x(1)), ind(x(2)), indt(t));
        %disp("recycled");
        %disp([ind(x(1)), ind(x(2))]);
    end
end

% Stretching, which is not necessary to consider in 2D
function dw = stretch(t, x)
    global w_t vx vy;
    %disp(ind(x(1)));
    %disp(ind(x(2)));
    h = 1;
    dw = abs(w_t(ind(x(1)), ind(x(2)), indt(t))/h) * [vx(ind(x(1))+1, ind(x(2))+1, indt(t)) - vx(ind(x(1))-1, ind(x(2))-1, indt(t)),
        vy(ind(x(1))+1, ind(x(2))+1, indt(t)) - vy(nd(x(1))-1, ind(x(2))-1, indt(t))
    ];
end

% Convert spatial coordinate to array index
function i = ind(x)
    global dx n_grid;
    i = round(x/dx) + 1;
    if i < 1
        i = 1;
    end
    if i > n_grid
    	i = n_grid;
    end
end

%Compute gradient of stream function using WoS
function grad = WoSgrad(t, x)
global WoS_eps n_grid dx;
	d = min([x(1), x(2), (n_grid - 1)*dx - x(1), (n_grid - 1)*dx - x(2)]);

    if d < WoS_eps
        grad = [0,0];
    else
        d_y = d*sqrt(rand);
        th = 2*pi*rand;
        th_y = 2*pi*rand;
        
        y = x + d_y*[cos(th_y), sin(th_y)];
        x_next = x + d*[cos(th), sin(th)];
        
        grad_g = (0.5/pi)*(y-x)*( 1/(norm(y-x)^2) + 1/(d^2));
        
        grad = (2/(d^2))*WoS(t, x_next)*(x_next-x) - vort(t, y)*grad_g;
        
    end
end


%Walk-on-Spheres method
function stream = WoS(t, x)
global WoS_eps n_grid dx;
	d = min([x(1), x(2), (n_grid - 1)*dx - x(1), (n_grid - 1)*dx - x(2)]);

    if d < WoS_eps
        stream = 0;
    else
        d_y = d*sqrt(rand);
        th = 2*pi*rand;
        th_y = 2*pi*rand;
        
        y = x + d_y*[cos(th_y), sin(th_y)];
        x_next = x + d*[cos(th), sin(th)];
        
        v = (4/3)*pi*(d_y^3);
        f = vort(t, [y(1), y(2)]);
        g = (1/(2*pi)) * log(d/norm(y-x));
        
        stream = WoS(t, x_next) + f*g*v;
        
    end
end

%Turn time value into index for array
function i = indt(t)
    global dt t_max;
    if t/dt > t_max
        i = round(t_max/dt) + 1;
        return;
    end
    i = round(t/dt)+ 1;
    if i < 1
        i = 1;
    end
end