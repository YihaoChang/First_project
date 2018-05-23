function plot_drag_lift(dt,nt,Cd,Cl,time_interval)

plot((0:nt)*dt,Cd,'linewidth',1.5);hold on;
plot((0:nt)*dt,Cl,'r--','linewidth',1.5);hold off;
set(gca,'fontsize',16);

interval = [time_interval,-0.5,1.5];
axis(interval)
% title('Drag and lift coefficient for stationary cylinder with Re = 100');
xlabel('Time')
ylabel('Drag and lift coefficient')
legend1=legend('Drag','Lift',2);
set(legend1,...
    'Position',[0.650446428571427 0.633630952380951 0.172321428571429 0.146726190476191]);
grid on

end