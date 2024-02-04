clear all;
close all;

% =========================================
% Partie I : Etude de la chaine DVBS-S    =
% dans un canal Gaussien                  =
%                                         =
%                                         =
% 1 Modulateur/ Démodulateur              =
%                                         =
% 1.1 Génération de l'enveloppe complexe  =
% associée au signal à transmettre        =
% &                                       =
% 1.2 Mise en place du récepteur          =
% en l'absence de canal                   =
% =========================================

% INITIALISATION
% Constantes
Fe = 24000;
Rb = 3000;
Tb = 1/Rb;


%% 1.1 Générarion de l'enveloppe complexe associée au signal à transmettre

% Modulation QPSK - 4 symboles
n = 2;    % nb de bits par symboles
M = 2^n;  % span 4 symbols car QPSK

Rs = Rb/n;
Ns = Fe/Rs; % nbe de symboles / échantillons
Ts = 1/Rs;
Te = 1/Fe;

% Génération aléatoire des n bits = information
nb_bits = 2000;
bits = randi([0 1], 1, nb_bits);

% Mapping des bits en symboles QPSK
% On passe les bits en +1 -1
bits =bits*2-1;

%  1 + j = 00 = -1 -1
% -1 + j = 10 =  1 -1
% -1 - j = 11 =  1  1
%  1 - j = 01 = -1  1

% Génération des symboles
symboles = zeros(1, nb_bits/2);
for i = 1:2:nb_bits
    symboles((i+1)/2) = -bits(i) - 1i * bits(1+i);
end



% Tracé de la constellation de l'émetteur
figure;
scatter(real(symboles), imag(symboles), 'filled');
title('Constellation QPSK (Avant transmission)');
xlabel('Partie réelle');
ylabel('Partie imaginaire');
grid on;

% Emetteur
signal = kron(symboles, [1 zeros(1, Ns -1)]);   % On prolonge de symbole avec N-1 zéros derrière

% Filtre de mise en forme 
% Racine de cosinus surélevé Beta = 0.35
Beta = 0.35;
h_n = rcosdesign(Beta, M, Ns,"sqrt"); 

% Signal en sortie de l'émetteur
x = conv(signal, h_n, 'full'); % Emission de notre signal 

% Tracé du signal modulé
figure;
plot(real(x), 'b');
hold on
plot(imag(x), 'r');
title('Signal modulé (partie réelle et imaginaire)');
xlabel('Temps');
ylabel('Amplitude');
legend('Partie réelle', 'Partie imaginaire');
grid on;
hold off;

% Pas de canal pour le moment



%% 1.2 Mise en place du récepteur en l'absence de canal
% Récepteur
% avec le filtre adapté
% donc on inverse la réponse impulsionnelle du filtre en émission
z = conv(x, fliplr(h_n), 'full');  % on inverse le filtre de mise en forme pour optenir le filtre adapté


% Tracé de la constellation en réception
figure;
scatter(real(z(128/2 + 1:Ns:end-128/2)), imag(z(128/2 + 1:Ns:end-128/2)), 'x');
title('Constellation QPSK (Après transmission)');
xlabel('Partie réelle');
ylabel('Partie imaginaire');
grid on;

% on enlève les 128 échantillons de retard
z = z(128/2 + 1:end - 128/2);

% Tracer du diagramme de l'oeil pour connaitre le moment optimal de la
% décision
eyediagram(z(Ns*2:end - Ns*2),Ns*2); % On regarde 2 symboles par fenetre

% Décision sur les parties réelles et imaginaires
% Reconstruction des symboles estimés à partir des décisions R et I
symboles_estimes = ((real(z(:,1:Ns:end)) > 0)*2-1) + 1i *((imag(z(:,1:Ns:end)) > 0)*2-1); 
% Décision sur la partie imaginaire et réel direct

% On se peut donc maitenant reconstituer l'information envoyé
bits_recu = zeros(1, nb_bits); % Initialisation des bits reçus

for i = 1:2:nb_bits
    bits_recu(i) = -real(symboles_estimes((i+1)/2));
    bits_recu(i+1) = -imag(symboles_estimes((i+1)/2));
end


% Calcul du taux d'erreur de symbole (TES)
mat_erreur_symb = symboles_estimes ~= symboles; 
nb_erreurs = sum(mat_erreur_symb); % Comptage des erreurs
taux_erreur_symbole = nb_erreurs / length(symboles); % Calcul du taux d'erreur

% Affichage du taux d'erreur
fprintf('TES = %f\n', taux_erreur_symbole);


% Calcul TEB 
mat_TEB = bits_recu ~= bits;
TEB = sum(mat_TEB)/nb_bits;
fprintf('TEB = %f\n', TEB);
