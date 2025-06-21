using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using System.Threading;
//using System.Threading.Tasks;
using System;
using Cysharp.Threading.Tasks;
using EyE.UnityAssetTypes;
using EyE.Threading;
namespace EyE.Geometry
{
    public class PolyhedronMaker : MonoBehaviour
    {
        public Mesh startingMesh;//optional: will be used during start to compute a starting polyhedron, IF startingPoly is null.
        public PolyhedronSO startingPoly;//optional: will be used as the starting polyhedron, if provided.

        public bool defaultToTetraNotIcosa;// if neither of the above two are selected- adefault platonic shape will be usedas the starting polyhedron

        public Polyhedron poly;
        public MeshFilter mf;
        public FacesAndNeighbors facesAndNeighborsOnMesh;
        public int numberOfIteration = 3;
        [Range(0, 0.5f)]
        public float truncateEdgeLengthFraction = 1f / 3f;
        public GameObject processingDisplayObject;

        private CommandLineProcessor commandLine;
        private Dictionary<string, Func<UniTask>> commandLineOperations;
        private void Awake()
        {
            commandLineOperations = new Dictionary<string, Func<UniTask>>()
                {
                    {"truncate", () => StartProcessingNow(() => poly.TruncAsync(truncateEdgeLengthFraction, cancelRef))},
                    {"dual", () => StartProcessingNow(() =>poly.DualAsync(cancelRef))},
                    {"tesselate", () => StartProcessingNow(() =>poly.TesselateTriangleByEdgeMiddlesAsync(cancelRef))},
                    {"radial", () => StartProcessingNow(() =>poly.TesselateFacesRadialAsync(cancelRef))},
                    {"spherize", () => StartProcessingNow(() =>poly.SpherizeAsync(1f, cancelRef))},
                    {"smooth", () => StartProcessingNow(() =>poly.SmoothAsync(cancelRef))}
                };

            commandLine = new CommandLineProcessor(commandLineOperations, cancelRef);
        }


        private void LaunchAsync<T>(Func<UniTask<T>> asyncFunc)
        {
#pragma warning disable CS4014
            StartProcessingNow(asyncFunc).Forget();
#pragma warning restore CS4014
        }

        private async UniTask StartProcessingNow<T>(Func<UniTask<T>> asyncFunc)
        {
            T result = await InternalLaunchAsync(asyncFunc);
            if (result is Polyhedron p) poly = p;
            else if (result is Mesh m)
            {
                await UniTask.SwitchToMainThread();
                SetDrawnMesh(m);
            }
        }

        private async UniTask<T> InternalLaunchAsync<T>(Func<UniTask<T>> asyncFunc)
        {
            await ProcessStart();
            T retVal = default(T);
            try
            {
                retVal = await asyncFunc();
            }
            catch (Exception ex)
            {
                Debug.LogError("Async error: " + ex);
                // Optionally: show error to user
            }
            finally
            {
                await ProcessDone();
            }
            return retVal;
        }
        bool isProcessing = false;
        async UniTask ProcessStart()
        {

            cancelRef.doCancel = false;
            await UniTask.SwitchToMainThread();
            isProcessing = true;
            await UniTask.SwitchToTaskPool();
        }
        async UniTask ProcessDone()
        {
            await UniTask.SwitchToMainThread();
            isProcessing = false;
        }

        async UniTask<Polyhedron> ResetPolyhedron()
        {
            await UniTask.SwitchToMainThread();
            facesAndNeighborsOnMesh = FacesAndNeighbors.CreateInstance<FacesAndNeighbors>();
            if (startingPoly != null)
            {
                poly = new Polyhedron(startingPoly.poly);
                SetDrawnMesh(startingPoly.mesh);
                facesAndNeighborsOnMesh = startingPoly.facesAndNeighbors;
            }
            else
            {
                if (startingMesh != null)
                {
                    SetDrawnMesh(startingMesh);
                    MeshData startingMeshData = new MeshData(startingMesh);
                    await UniTask.SwitchToThreadPool();
                    Polyhedron p = await Polyhedron.CreateFromMesh(startingMeshData, cancelRef);
                    await p.ToMeshAsync(facesAndNeighborsOnMesh);//output will be functionally identical to starting mesh, except perhaps uvs
                    //generate mesh output currently ignored, in favor of startingmesh, but not sure if it should...
                }
                else
                {
                    if (defaultToTetraNotIcosa)
                        poly = Polyhedron.ComputeTetrahedron();
                    else
                        poly = Polyhedron.ComputeIcosahedron();
                    SetDrawnMesh(poly.ToMesh(facesAndNeighborsOnMesh));
                }
            }
            return poly;
        }

        // Start is called before the first frame update
        void Start()
        {
            LaunchAsync(() => ResetPolyhedron());
        }

        void Update()
        {
            processingDisplayObject.SetActive(isProcessing);
        }

        CancelBoolRef cancelRef = new CancelBoolRef();
        // Method to cancel the task
        void CancelTask()
        {
            cancelRef.doCancel = true;
            cancelRef = new CancelBoolRef();
        }

        public Gradient orderGradient = new Gradient();


        void SetDrawnMesh(Mesh mesh)
        {
            mf.sharedMesh = mesh;
            MeshCollider mc;
            if (TryGetComponent<MeshCollider>(out mc))
            {

                mc.sharedMesh = mesh;
            }
        }

        public void SaveSOs()
        {
            string namePrefix = poly.faces.Count.ToString();
            if (poly.faces[0].corners.Count == 3)
                namePrefix += "Tri";
            Mesh meshRef = mf.sharedMesh;
            if (meshRef == null)
            {
                Debug.LogWarning("Unable to save mesh asset: unable to find one assigned to meshFilter.sharedMesh");
            }
            else
            {
                string filename = "assets/" + namePrefix + "Mesh.asset";
                if (UnityEditor.AssetDatabase.LoadAssetAtPath<UnityEngine.Object>(filename) != null)
                    UnityEditor.AssetDatabase.DeleteAsset(filename);
                UnityEditor.AssetDatabase.CreateAsset(meshRef, filename);
            }
            FacesAndNeighbors savedInstance = FacesAndNeighbors.Instantiate<FacesAndNeighbors>(facesAndNeighborsOnMesh);
            savedInstance.meshRef = mf.sharedMesh;
            string facesFilename = "Assets/" + namePrefix + "FacesAndNeighbors.asset";
            if (UnityEditor.AssetDatabase.LoadAssetAtPath<UnityEngine.Object>(facesFilename) != null)
                UnityEditor.AssetDatabase.DeleteAsset(facesFilename);
            UnityEditor.AssetDatabase.CreateAsset(savedInstance, facesFilename);

            PolyhedronSO polyhedronSO = PolyhedronSO.CreateInstance<PolyhedronSO>();
            polyhedronSO.poly = poly;
            polyhedronSO.facesAndNeighbors = savedInstance;
            polyhedronSO.mesh = mf.sharedMesh;
            string polyFilename = "Assets/" + namePrefix + "Polyhedron.asset";
            if (UnityEditor.AssetDatabase.LoadAssetAtPath<UnityEngine.Object>(polyFilename) != null)
                UnityEditor.AssetDatabase.DeleteAsset(polyFilename);
            UnityEditor.AssetDatabase.CreateAsset(polyhedronSO, polyFilename);

        }
        private string commandInput = "";
        private bool recordMode = false;

        void OnGUI()
        {
            GUI.enabled = !isProcessing;

            GUILayout.Label("Command Line:");
            commandInput = GUILayout.TextField(commandInput, GUILayout.Width(400));

            if (GUILayout.Button(recordMode ? "Stop Recording" : "Record Commands"))
            {
                recordMode = !recordMode;
            }

            GUI.enabled = !string.IsNullOrWhiteSpace(commandInput);
            if (GUILayout.Button("Execute Command Line"))
            {
                commandLine.ExecuteSequence(commandInput)
                    .ContinueWith(() => StartProcessingNow(()=>poly.ToMeshAsync(facesAndNeighborsOnMesh)))
                    .Forget();//fire and forget sequence
            }

            GUI.enabled = !isProcessing;
            GUILayout.Space(20);
            GUILayout.Label("Operations:");

            foreach (var kvp in commandLineOperations)
            {
                string cmd = kvp.Key;
                if (GUILayout.Button(cmd))
                {
                    if (recordMode)
                        commandInput += cmd + " ";
                    else
                        kvp.Value.Invoke();
                }
            }

            GUI.enabled = isProcessing;
            if (GUILayout.Button("Cancel"))
            {
                CancelTask();
            }
            GUI.enabled = true;
            if (GUILayout.Button("Reset"))
            {
                LaunchAsync(() => ResetPolyhedron());
            }
            if (GUILayout.Button("IcoShpere"))
            {
                (Polyhedron, Mesh) ico = poly.Icosphere(numberOfIteration);
                poly = ico.Item1;
                SetDrawnMesh(ico.Item2);
            }
            if (GUILayout.Button("RecomputePolyMesh"))
            {
                facesAndNeighborsOnMesh = FacesAndNeighbors.CreateInstance<FacesAndNeighbors>();
                LaunchAsync<Mesh>(() => poly.ToMeshAsync(facesAndNeighborsOnMesh));
            }
            if (GUILayout.Button("CheckPoly"))
            {
                Debug.Log("Checking Polyhedron Integrity/Regularity");
                poly.CheckIntegrity(transform, orderGradient);
            }
            if (GUILayout.Button("Save"))
            {
                SaveSOs();
            }

            if (isProcessing)
            {
                DrawOscillatingBox(new Vector2(Screen.width / 2, 500));
            }

            void DrawOscillatingBox(Vector2 pos)
            {
                float frac = Mathf.Sin(Time.time * speed) * 0.4f + 0.6f;
                float width = boxSize * (1f + frac);
                float height = boxSize * 0.5f;

                Rect boxRect = new Rect(pos.x - (width / 2f), pos.y, width, height);
                GUIStyle noBackgroundStyle = new GUIStyle(GUI.skin.textArea);
                noBackgroundStyle.normal.background = Texture2D.blackTexture;
                noBackgroundStyle.border = new RectOffset(0, 0, 0, 0);
                noBackgroundStyle.normal.textColor = Color.white;
                noBackgroundStyle.padding = new RectOffset(0, 0, 0, 0);
                noBackgroundStyle.alignment = TextAnchor.MiddleCenter;

                GUI.Box(boxRect, "");
                GUI.TextArea(new Rect(pos.x - (boxSize / 2f), pos.y - height / 4f, 200, 60), "Processing", noBackgroundStyle);
            }
        }

        /*
        void OnGUI()
        {
            GUI.enabled = !isProcessing;
            *//*
            GUILayout.Label($"Number of Iterations: {numberOfIteration}");
            numberOfIteration = Mathf.RoundToInt(GUILayout.HorizontalSlider(numberOfIteration, 1, 10));

            if (GUILayout.Button(new GUIContent("Multiple TriangleTesselation, Spherize, Dual Once", "Use on polyhedrons with triangular faces")))
                StartCoroutine(DoMultipleTriEdgeTessSpherizeDual());
            if (GUILayout.Button(new GUIContent("Multiple RadialFaceTesselation, Spherize, Dual Once", "Use on polyhedrons with triangular faces")))
                StartCoroutine(DoRadialFaceTessSpherizeDual());


            if (GUILayout.Button(new GUIContent("Multiple TriangleTesselation then, Spherize Once", "Use on polyhedrons with triangular faces")))
                StartCoroutine(DoTessThenSpherizeOnce());
            *//*
            if (GUILayout.Button("TruncatePoly"))
            {
                LaunchAsync(() => poly.TruncAsync(truncateEdgeLengthFraction, cancelRef));
            }
            if (GUILayout.Button("DualPoly"))
            {
               LaunchAsync(() => poly.DualAsync(cancelRef));
            }
            if (GUILayout.Button("Triangle Edge Tesselate"))
            {
                LaunchAsync(() => poly.TesselateTriangleByEdgeMiddlesAsync(cancelRef));
            }
            if (GUILayout.Button("Radial Face Tesselate"))
            {
                LaunchAsync(() => poly.TesselateFacesRadialAsync(cancelRef));
            }

            if (GUILayout.Button("Spherize"))
            {
                LaunchAsync(() => poly.SpherizeAsync(1f, cancelRef));
            }
            if (GUILayout.Button("Smooth"))
            {
                LaunchAsync(() => poly.SmoothAsync(cancelRef));
            }
            if (GUILayout.Button("IcoShpere"))
            {
                (Polyhedron, Mesh) ico = poly.Icosphere(numberOfIteration);
                poly = ico.Item1;
                SetDrawnMesh(ico.Item2);
            }
            if (GUILayout.Button("RecomputePolyMesh"))
            {
                facesAndNeighborsOnMesh = FacesAndNeighbors.CreateInstance<FacesAndNeighbors>();
                LaunchAsync<Mesh>(() => poly.ToMeshAsync(facesAndNeighborsOnMesh));
            }
            if (GUILayout.Button("CheckPoly"))
            {
                Debug.Log("Checking poly");
                poly.CheckIntegrity(transform, orderGradient);
            }
            if (GUILayout.Button("Save"))
            {
                SaveSOs();
            }


            // Reset button to reset startingPoly
            if (GUILayout.Button("Reset"))
            {
                LaunchAsync(() => ResetPolyhedron());
            }
            GUI.enabled = isProcessing;
            if (GUILayout.Button("Cancel"))
            {
                CancelTask();
            }
            GUI.enabled = true;


            if (isProcessing)
            {

                DrawOscillatingBox(new Vector2(Screen.width / 2, 500));
            }


            void DrawOscillatingBox(Vector2 pos)
            {
                // Update box size using a sine wave
                float frac = Mathf.Sin(Time.time * speed) * 0.4f + 0.6f; // Oscillates between 0.2 and 1

                // Calculate dimensions
                float width = boxSize * (1f + frac);
                float height = boxSize * 0.5f;// * (1f - frac); // Opposite phase

                // Draw the box centered on the screen
                Rect boxRect = new Rect(pos.x - (width / 2f), pos.y, width, height);
                GUIStyle noBackgroundStyle = new GUIStyle(GUI.skin.textArea);
                noBackgroundStyle.normal.background = Texture2D.blackTexture;// null;
                noBackgroundStyle.border = new RectOffset(0, 0, 0, 0);
                noBackgroundStyle.normal.textColor = Color.white; // or any color you want
                noBackgroundStyle.padding = new RectOffset(0, 0, 0, 0);
                noBackgroundStyle.alignment = TextAnchor.MiddleCenter;
                GUI.Box(boxRect, "");//, "Processing...");
                GUI.TextArea(new Rect(pos.x - (boxSize / 2f), pos.y - height / 4f, 200, 100), "Processing", noBackgroundStyle);
            }
        }
        */
        public float boxSize = 200f; // Box size multiplier (oscillates)
        public float speed = 1f;     // Oscillation speed
    }
}