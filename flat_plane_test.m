%Saves rerunning model

query = questdlg('run model?');
switch query
    case 'Yes'   
    run("scatteredPressure_main.m")
    case 'No'
end

lvals = [];
mphopen('flat_plane')
figure()
hold on
for ii = 1:2:(Probe.domain - 1)
    Probe.radius = ii; %radius of arc in meters
    Probe.Coordinates(1,:) = Probe.radius*cos(Probe.theta_vector); %Probe x coordinates
    Probe.Coordinates(2,:) = Probe.radius*sin(Probe.theta_vector); %Probe y coordinates
    
    Psflatnum = mphinterp(model,'acpr.p_s','coord',Probe.Coordinates);
    
    n_d = length(Probe.theta_vector);

    % Flat plane
    SIflatnum = abs(Psflatnum).^2; %sound intensity
    SIsumflatnum= sum(SIflatnum,2);
    SIsqflatnum = sum(SIflatnum.^2,2);
    delta_flatnum = (SIsumflatnum.^2 - SIsqflatnum)./((n_d-1)*(SIsqflatnum));
    
    name = "r =" + string(Probe.radius) + " m";
    plot(Freq.Vector,delta_flatnum,"LineWidth",1,"DisplayName",name)
   
end  

plot(Freq.Vector,deltaf,"LineWidth",1,"LineStyle","--","DisplayName","TMM")
lgd = legend;
lgd.NumColumns = 2;
legend('Location','southeast')
title(['parametric study of effect of probe distance - domain radius: ',string(Probe.domain),' m'])
ylim([0.85, 1])
xlabel("Hz")
ylabel("diffusion coefficient")