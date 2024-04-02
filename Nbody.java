public class Nbody {
    public static void main(String[] args) {
        /*Accepting total time and length of the increment. This will come into play when we calculate forces*/
        double T= Double.parseDouble(args[0]);
        double deltaT = Double.parseDouble(args[1]);
        /*Gravitational constant*/
        double G = 6.67E-11;
        /*Next we move onto collecting the inputs. First, since they are single tokens, we read the number of bodies and the universe
        radius into standard input separately before collecting individual body details. */
        int n = StdIn.readInt();
        double radius = StdIn.readDouble();
        /*Reading all of the details of all the bodies into a single array, bodyDetails, of strings.*/
        /**@param bodyDetails is an array consisting of all the information in the input file's table read in as strings*/
         /** @param numerics is an n by 5 array that converts the string array bodyDetails(which are actually number)*/
         /*into a two-dimensional array of doubles, which we can now usefully parse for future use*/
        /*(6*i)+j is a command we use to skip over the picture files for each body. This works because they occur every sixth element */
        String[] bodyDetails = StdIn.readAllStrings();
        Double[][] numerics = new Double[n][5];
        for (int i = 0; i <= n - 1; i++) {
            int j = 0;
            while (j <= 4) {
                numerics[i][j] = Double.parseDouble(bodyDetails[(6 * i) + j]);
                j += 1;
            }
        }
        /**@param mass is a 1-dimensional array that parses through all the elements of numerics and assembles all the body masses into one array* */
        Double[] mass = new Double[n];
        for (int k = 0; k <= n - 1; k++) {
            mass[k] = numerics[k][4];
        }
        /*System.arraycopy allows one to effectively assemble new arrays from an existing one by slicing it and assigning the sliced section
        * to another array.We assemble arrays containing the x, y coordinates of position and velocity of each body in this way  */
        Double[][] position = new Double[n][2];
        Double[][] velocity = new Double[n][2];
        for (int i = 0; i <= n - 1; i++) {
            System.arraycopy(numerics[i], 0, position[i], 0, 2);
            System.arraycopy(numerics[i], 2, velocity[i], 0, 2);
        }
        Double[][] force = new Double[n][2];
        /*Now that we have created sub-arrays with position and velocity components of each body, we can begin the computation*/
        /*Setting the drawing parameters, like the background and scale. as well as setting the time to start from 0.0.
        * We enable DoubleBuffering to draw all our objects at once for each iteration of the overall loop,
        * The next three loops set all forces at the current position to 0.0, calculate the x-components and calculate the y components
        * using Newton's formula respectively.
        * Next we update each velocity in x and y directions by adding to each existing term,
        * the product of its acceleration with the time interval. Acceleartion is included in our calculation by dividing by mass each time
        * Accordingly we update each x, y coordinate of the position.
        * Each body is then drawn at this new set of positions and then the frame is displayed with the show command*/

        StdDraw.setScale(-radius,radius);
        double elapsedT = 0.0;
        while (elapsedT <= T) {
            StdDraw.enableDoubleBuffering();
            StdDraw.picture(0,0,"starfield.jpg");
            for(int b=0; b<n; b++){
                for(int c=0;c<=1;c++){
                    force[b][c]=0.0;
                }
            }
            for(int k=0;k<n; k++){
                for(int l=0;l<n;l++){
                    if(l!=k){
                        force[k][0]+=(G*mass[k]*mass[l]/Math.pow(Math.pow(position[k][0]-position[l][0],2)+Math.pow(position[k][1]-position[l][1],2),1.5))*(position[l][0]-position[k][0]);
                    }
                }
            }
            for(int k=0;k<n; k++){
                for(int l=0;l<n;l++){
                    if(l!=k){
                        force[k][1]+=(G*mass[k]*mass[l]/Math.pow(Math.pow(position[k][0]-position[l][0],2)+Math.pow(position[k][1]-position[l][1],2),1.5))*(position[l][1]-position[k][1]);
                    }
                }
            }
            for(int k=0; k<n;k++){
                velocity[k][0] += ((force[k][0]/mass[k]) * deltaT);
            }
            for(int k=0; k<n;k++){
                velocity[k][1] += ((force[k][1]/mass[k]) * deltaT);
            }
            for(int k=0;k<n;k++){
                position[k][0] += velocity[k][0]*deltaT;
            }
            for(int k=0; k<n;k++){
                position[k][1] += velocity[k][1]*deltaT;
            }
            for (int k= 0; k<n; k++) {
                StdDraw.picture(position[k][0], position[k][1], bodyDetails[(6 * k) + 5]);
            }
            StdDraw.show();
            StdDraw.clear();
            elapsedT+=deltaT;
           }
        System.out.println(n);
        System.out.println(radius);
        for(int i=0;i<n;i++){
            System.out.println(position[i][0]+" "+position[i][1]+" "+velocity[i][0]+" "+velocity[i][1]+" "+bodyDetails[6*i+4]+" "+bodyDetails[6*i+5]);
        }
       }
   }


