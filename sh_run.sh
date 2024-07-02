# Parameters
E0=1.0
omega0=1.0
Tprop=600
dt=0.01
omega_c=0.5
eta=0.1
rx=1.0
rz=0.0

# Born (eq.)
echo 'Born (eq.)'
echo "'Born'" > inp_born_eq
echo '0d0' >> inp_born_eq
echo $omega0 >> inp_born_eq
echo $Tprop >> inp_born_eq
echo $dt >> inp_born_eq
echo $omega_c >> inp_born_eq
echo $eta >> inp_born_eq
echo $rx >> inp_born_eq
echo $rz >> inp_born_eq

cp inp_born_eq inp
./test < inp > log_born_eq.log

cp pop_t_born.out pop_t_born_eq.out


# Born
echo 'Born'
echo "'Born'" > inp_born
echo $E0 >> inp_born
echo $omega0 >> inp_born
echo $Tprop >> inp_born
echo $dt >> inp_born
echo $omega_c >> inp_born
echo $eta >> inp_born
echo $rx >> inp_born
echo $rz >> inp_born

cp inp_born inp
./test < inp > log_born.log


# Redfield (eq.)
echo 'Redfield (eq.)'
echo "'Redfield'" > inp_redfield_eq
echo '0d0' >> inp_redfield_eq
echo $omega0 >> inp_redfield_eq
echo $Tprop >> inp_redfield_eq
echo $dt >> inp_redfield_eq
echo $omega_c >> inp_redfield_eq
echo $eta >> inp_redfield_eq
echo $rx >> inp_redfield_eq
echo $rz >> inp_redfield_eq

cp inp_redfield_eq inp
./test < inp > log_redfield_eq.log

cp pop_t_redfield.out pop_t_redfield_eq.out


# Redfield
echo 'Redfield'
echo "'Redfield'" > inp_redfield
echo $E0 >> inp_redfield
echo $omega0 >> inp_redfield
echo $Tprop >> inp_redfield
echo $dt >> inp_redfield
echo $omega_c >> inp_redfield
echo $eta >> inp_redfield
echo $rx >> inp_redfield
echo $rz >> inp_redfield

cp inp_redfield inp
./test < inp > log_redfield.log



# Lindblad (eq.)
echo 'Lindblad (eq.)'
echo "'Lindblad'" > inp_lindblad_eq
echo '0d0' >> inp_lindblad_eq
echo $omega0 >> inp_lindblad_eq
echo $Tprop >> inp_lindblad_eq
echo $dt >> inp_lindblad_eq
echo $omega_c >> inp_lindblad_eq
echo $eta >> inp_lindblad_eq
echo $rx >> inp_lindblad_eq
echo $rz >> inp_lindblad_eq

cp inp_lindblad_eq inp
./test < inp > log_lindblad_eq.log

cp pop_t_lindblad.out pop_t_lindblad_eq.out

# Lindblad
echo 'Lindblad'
echo "'Lindblad'" > inp_lindblad
echo $E0 >> inp_lindblad
echo $omega0 >> inp_lindblad
echo $Tprop >> inp_lindblad
echo $dt >> inp_lindblad
echo $omega_c >> inp_lindblad
echo $eta >> inp_lindblad
echo $rx >> inp_lindblad
echo $rz >> inp_lindblad

cp inp_lindblad inp
./test < inp > log_lindblad.log


# Analysis
echo 'Analysis'
echo "'analysis'" > inp_analysis
echo $E0 >> inp_analysis
echo $omega0 >> inp_analysis
echo $Tprop >> inp_analysis
echo $dt >> inp_analysis
echo $omega_c >> inp_analysis
echo $eta >> inp_analysis
echo $rx >> inp_analysis
echo $rz >> inp_analysis

cp inp_analysis inp
./test < inp > log_analysis.log
