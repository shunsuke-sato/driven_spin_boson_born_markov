# Born
echo 'Born'
cp inp_born_eq inp
./test < inp > log_born_eq.log

cp pop_t_born.out pop_t_born_eq.out

cp inp_born inp
./test < inp > log_born.log

# Redfield
echo 'Redfield'
cp inp_redfield_eq inp
./test < inp > log_redfield_eq.log

cp pop_t_redfield.out pop_t_redfield_eq.out

cp inp_redfield inp
./test < inp > log_redfield.log

# Lindblad
echo 'Lindblad'
cp inp_lindblad_eq inp
./test < inp > log_lindblad_eq.log

cp pop_t_lindblad.out pop_t_lindblad_eq.out

cp inp_lindblad inp
./test < inp > log_lindblad.log


# Analysis
echo 'Analysis'
cp inp_analysis inp
./test < inp > log_analysis.log
