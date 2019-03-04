% GMM
load sim_data.mat;
N = 75000; iter = 1;
U = rand(N,2); param_mat = zeros(7,iter);
theta0 = [10;-4;-0.8;-1.5;-3;-1;0.8];
options = optimset('Display', 'iter', 'PlotFcns', @optimplotx, 'MaxIter', 1200, 'OutputFcn', @outfun);
for k = 1:iter
    k
    gen_data;
    param_mat(:,k) = fminsearch(@(param)mom_dist(param, est_data, U, z2zz), theta0, options);
end

function stop = outfun(x, optimValues, state)
stop = false;
hold on;
dlmwrite('x.txt',x);
end

function d = mom_dist(param, data, U, z2zz)
    zone = [1:5,7:16]; N = size(data,1); choice = data(:,1);
    true_mkt_share_taxi = zeros(15,24); zz = z2zz(data(:,3),2);
    true_mkt_share_uber = zeros(15,24);
    h_data = floor((data(:,2)-1)/6);
    for z_id = 1:15
        for h = 0:23
            true_mkt_share_taxi(z_id,h+1) = sum(choice==1 & zz==zone(z_id) & h_data==h)/sum(zz==zone(z_id) & h_data==h);
            true_mkr_share_uber(z_id,h+1) = sum(choice==2 & zz==zone(z_id) & h_data==h)/sum(zz==zone(z_id) & h_data==h);
        end
    end
    true_mom = [reshape(true_mkt_share_taxi,360,1); reshape(true_mkr_share_uber,360,1)];
    
    alpha = param(1:2); beta = param(3:4); gamma = param(5:6); sigma = param(7);
    z1 = data(:,3); z2 = data(:,4); taxi_P = data(:,7); taxi_T = data(:,8); uber_P = data(:,9); 
    dist = data(:,5); uber_T = data(:,10); mta_P = data(:,11); mta_T = data(:,12);
    
    U_taxi = beta(1)*taxi_P + gamma(1)*taxi_T;
    U_uber = beta(1)*uber_P + gamma(1)*uber_T;
    U_nest = [U_taxi,U_uber]-max([U_taxi,U_uber],2);
    U_ride = alpha(1) + alpha(2) * (dist > 3) + max([U_taxi,U_uber]) + sigma(1)*log(sum(exp(U_nest),2));
    ride_prob = exp(U_nest)./sum(exp(U_nest),2);
    ride_choice = sum(bsxfun(@gt, U(:,1), cumsum(ride_prob,2)),2)+1;
    U_mta =  beta(2)*mta_P + gamma(2)*mta_T;
    U_vect = [U_ride,U_mta]-max([U_ride,U_mta],2);
    nest_prob = exp(U_vect)./sum(exp(U_vect),2);
    nest_choice = sum(bsxfun(@gt, U(:,2), cumsum(nest_prob,2)),2)+1;
    sim_choice = zeros(N,1);
    sim_choice(nest_choice==2) = 3;
    sim_choice(nest_choice==1) = ride_choice(nest_choice==1);
    
    sim_mkt_share_taxi = zeros(15,24);
    sim_mkt_share_uber = zeros(15,24);
    for z_id = 1:15
        for h = 0:23
            sim_mkt_share_taxi(z_id,h+1) = sum(sim_choice==1 & zz==zone(z_id) & h_data==h)/sum(zz==zone(z_id) & h_data==h);
            sim_mkr_share_uber(z_id,h+1) = sum(sim_choice==2 & zz==zone(z_id) & h_data==h)/sum(zz==zone(z_id) & h_data==h);
        end
    end
    sim_mom = [reshape(sim_mkt_share_taxi,360,1); reshape(sim_mkr_share_uber,360,1)];
    d = norm(sim_mom-true_mom);
end
