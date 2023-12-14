
name="sub_theta"
sh_command="gfortran -shared "${name}".f90 -o "${name}".so"

echo ""
$sh_command

if [ $? -eq 0 ]; then
    echo "Library created successfully"
else
    echo "Error in compiling"
fi

